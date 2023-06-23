%% FDA library drug screen in e. coli
%data analysis and plotting
% May 21, 2021
% last updated June 06, 2021
% 16 plate screen of del-ybcN (~WT e. coli) against 1280 drugs
% By Mariana & Brittany 
%% User definitions
filetosave = 'megadrugscreen3.mat'; %% name of the matlab file storing 
% analysis variables created in this script
xfiletosave = '20210603sortedAntihostDrugs.xlsx';%% name of the excel file 
% holding drugs metadata, normalized OD sorted by increasing normalized OD values
saveFigs = false; % if true save figures created in analysis

flagTerm = {'antibacterial','antbacterial','antimicrobial','antibiotic','antiseptic','antiinfect'}; % terms to flag as antibiotics
cutoff_STD_multiplier = 1; % hits are detected by comparing their OD to the 
% OD distribution of the control, This multiplier is for how many STD units to use

nConditions = 3; % types of drugs in the screen: DMSO, antibacterial, antihost
cont = 1; % index of control drugs in cells
ab = 2; % index of antibiotic drugs in cells
ah = 3; % index of antihost drugs in cells
legendConditions = {'control';'antibiotic';'antihost'};
colorConditions = [0 1 0; 1 0 0; 0 0 1];

% choose number of timepoints to compare hits
inx_tp = [2 3]; % index of time points of interest
num_tp = 2; % number of time points of interest
label_tp = {'9 hours','11 hours'}; % labels of time points of interest
% label_tp = {strcat(num2str(inx_tp(1)-1),'hour'), strcat(num2str(inx_tp(2)-1),'hour'),strcat(num2str(inx_tp(3)-1),'hour'),};

% index of the timepoint of OD we chose as optimal
ODnormTime = 2;
% ODnormTime = 10;
%% load in data  - parced and normalized in seperate file
load drugscreen_normOD.mat % the name of Screen 1 
%from original source
% Definition of loaded variables %
% Drugs = drug metadata
% normODzero = background removed OD data
% normOD = normalized OD data (divided by median DMSO OD per plate)
% norm_OD_plot = combined table of drug metadata and normalized OD reads
% time = timepoints
% idx_DMSO = boolean index of DMSO well positions
%% Labeling antimicrobial drugs
% makes true false array for antimicrobial drugs, used to index
% different drug types on the next section
annotationArray = Drugs(:,end); % Drug annotation
drugArray = Drugs(:,2); % drug Names
nDrugs = length(Drugs);
% analyze the drug annotation to identify antimicrobial drugs 
tf = false(size(annotationArray)); % a true-false array to mark antimicrobial drugs
for i=1:length(annotationArray)
    curAnnotation = annotationArray{i};
    if(~isnan(curAnnotation)) % testing for flagTerm only if the current annotation string is not empty
        for j=1:length(flagTerm) % iterate over all flagTerm to find a match
            if(~isempty(strfind(curAnnotation,flagTerm{j}))) % true if the term is found somewhere in curAnnotation 
                tf(i) = true;
            end
        end
    end
end

%% Creates indices for each defined condition
inx_all = [1:nDrugs]'; 
inx_DMSO = find(strcmp(drugs, 'DMSO')); % Index of DMSO control well positions (boolean vector)
% store indexes in a single matrix; column1: control, column2:
% antimicrobial, colummn3: antihost
inx_conditions = {}; % create a cell to store indexes of all drug types
inx_conditions{1,cont} = setdiff(inx_DMSO,[1:96:nDrugs]); % Index of DMSO control well positions (boolean vector) This also removes A1

inx_drugs = inx_all(~ismember(inx_all,inx_DMSO)); % any well that is not a control well
inx_conditions{1,ab} = find(tf); % any well that is flagged as antimicrobial by the flag terms
inx_conditions{1,ah} = inx_drugs(~ismember(inx_drugs,inx_conditions{1,ab})); % any drug well that is not antibiotics

clear inx_all

%% Define a cutoff by the OD distribution of the control wells
cutoffOD = zeros(size(normOD,2),1); % empty array to hold the cutoff values
for i=1:size(normOD,2) %% iterate over all time points 
    curContOD = normOD(inx_conditions{1,cont},i); % get the normOD of all control wells
    curContOD_med = median(curContOD); 
    curContOD_std = std(curContOD);
    cutoffOD(i) = curContOD_med-cutoff_STD_multiplier*curContOD_std; % cutoff is defined as X-std away from the media OD on control wells
end

%% Find hits for each time point of the screen
nTime = size(normOD,2); % number of time points
% create a cell containing a true-false array: is the drug a hit by each
% time point, for each drug type
tf_hits{1,cont} = false(nDrugs,nTime);
tf_hits{1,ab} = false(nDrugs,nTime);
tf_hits{1,ah} = false(nDrugs,nTime);

n_hits = zeros(nTime,nConditions); % creates a matrix with the number of 
% hits for each time point (rows), for each drug type (columns)

% iterate over all time points and find all hits, store them in cell arrays
for iCond = 1:nConditions
    for  i=1:nTime
        tf = normOD(inx_conditions{1,iCond},i)<cutoffOD(i);
        n_hits(i,iCond) = sum(tf);
        tf_hits{1,iCond}(inx_conditions{1,iCond}(tf),i) = true;
    end
end

%% plot the histogram of normOD by time point
h = figure('Color','white'); hold on;
edges = [-0.11:0.02:1.51];
inx_tp = [2,3]; % time points to inspect
n = length(inx_tp); % useful variable for later

for iCond = 1:nConditions
    for i=1:n
        cur_inx_tp = inx_tp(i); % the the index of the time point for the current iteration
        subplot(3,n,n*(iCond-1)+i); hold on;
        histogram(normOD(inx_conditions{1,iCond},cur_inx_tp),edges,'normalization','probability','facecolor',colorConditions(iCond,:));
        xline(cutoffOD(cur_inx_tp),'k',{'-1 Ctrl Standard dev'},'LineWidth',1); % plot a line to show the cutoff
        n_drugs = length(inx_conditions{1,iCond}); % the number of wells in the distribution
        cur_n_hits = sum(tf_hits{1,iCond}(:,cur_inx_tp)); % number of toxic drug by the cutoff
        strTitle = strcat(legendConditions{iCond,1},'(', num2str(cur_n_hits), '/', num2str(n_drugs), ')'); % we build a string with the values we just calculated
        title(strTitle); grid on; box on;
    end
end
set(h,'Position',[551,273,537,348]); % adjust the figure size for 2
% timepoints
% set(h,'Position',[28 269 1364 486]);% adjust the figure size for 12
% timepoints
if saveFigs
    saveas(h,'normalized_OD_by_drug_type_decreased_growth_hits', 'png');
end
%% plot the histogram of normOD by time point for the paper - ADDED 01/20/2023 by MN
h = figure('Color','white'); hold on;
edges = [-0.11:0.02:1.51];
inx_tp = [3]; % time points to inspect
n = length(inx_tp); % useful variable for later
subplot(2,1,1);
histogram(normOD(inx_conditions{1,cont},inx_tp),edges,'normalization','probability','facecolor',[0.8 0.8 0.8]);
xline(cutoffOD(cur_inx_tp),'k',{'-1 Ctrl Standard dev'},'LineWidth',1); % plot a line to show the cutoff
n_drugs = length(inx_conditions{1,cont}); % the number of wells in the distribution
cur_n_hits = sum(tf_hits{1,cont}(:,inx_tp)); % number of toxic drug by the cutoff
strTitle = strcat(legendConditions{cont,1},'(', num2str(cur_n_hits), '/', num2str(n_drugs), ')'); % we build a string with the values we just calculated
title(strTitle); grid on; box on;
subplot(2,1,2);
histogram(normOD([inx_conditions{1,ab}; inx_conditions{1,ah}],inx_tp),edges,'normalization','probability','facecolor',[0 0 0]);
xline(cutoffOD(cur_inx_tp),'k',{'-1 Ctrl Standard dev'},'LineWidth',1); % plot a line to show the cutoff
n_drugs = numel([inx_conditions{1,ab}; inx_conditions{1,ah}]); % the number of wells in the distribution
cur_n_hits = sum([tf_hits{1,ab}(:,inx_tp); tf_hits{1,ah}(:,inx_tp)]); % number of toxic drug by the cutoff
strTitle = strcat('ab + ah (', num2str(cur_n_hits), '/', num2str(n_drugs), ')'); % we build a string with the values we just calculated
title(strTitle); grid on; box on; set(gcf,"Position",[476 446 560 420]);
saveas(h,[date 'normalized_OD_Fig1B.png']);
%%
h = figure('Color','white'); hold on;
edges = [-0.11:0.02:1.51];
inx_tp = [3]; % time points to inspect
n = length(inx_tp); % useful variable for later
histogram(normOD(inx_conditions{1,cont},inx_tp),edges,'normalization','probability','facecolor',[0.8 0.8 0.8]);
xline(cutoffOD(cur_inx_tp),'k',{'-1 Ctrl Standard dev'},'LineWidth',1); % plot a line to show the cutoff
histogram(normOD([inx_conditions{1,ab}; inx_conditions{1,ah}],inx_tp),edges,'normalization','probability','facecolor',[0 0 0]);
grid on; box on;
%%%%%%%%%%%%%%%%%%%%%%%
%% Creates a list of drug hits
hitNames = {};
hitAnnotations = {};
for iCond = 1:nConditions
    for i=1:nTime;
        cur_tf_hits = tf_hits{1,iCond}; % extracts the hits for condition i
        cur_inx_hits = find(cur_tf_hits(:,i)); % finds the position of each hit in condition i
        hitNames{1,iCond}{1,i} = drugArray(cur_inx_hits); % this variable will hold the hit names
        hitAnnotations{i,iCond} = annotationArray(cur_inx_hits); % this variable will hold the hit annotations
    end
end

%% plot a Venn diagram for overlap in hits
comparison = nchoosek(inx_tp,2); % coordinates for conditions to compare
nHits = zeros(numel(inx_tp),nConditions-1); % will store total number of 
% hits at each time point
intHits = zeros(length(comparison)+1,nConditions-1); % will store total 
% number of shared hits for each comparison
for iCond=2:nConditions
    cur_tf_hits = tf_hits{1,iCond};
    nHits(:,iCond-1) = sum(cur_tf_hits(:,inx_tp))'; % counts the number of 
    % hits in each condtion(column) at desired time points (row)
    for iCompare=1:size(comparison,1)
        intHits(iCompare,iCond-1) = sum(sum([cur_tf_hits(:,comparison(iCompare,:))],2) == 2);
        % common hits in combinations of time points specified by comparison
    end
    intHits(end,iCond-1) = sum(sum([cur_tf_hits(:,inx_tp)],2) == num_tp); % finds
    % common hits across all time points
end

for iCond = 2:nConditions
    h1 = figure; hold on;
    venn(nHits(:,iCond-1),intHits(:,iCond-1));
    legend(label_tp);
    if num_tp == 2
        tleg = strcat(label_tp{1},'-', label_tp{2},'-',legendConditions{iCond}, '-hits');
    elseif num_tp == 3
        tleg = strcat(lab{1},'-', lab{2},'-', lab{3},'-',legendConditions{iCond}, '-hits');
    end
    title(tleg)
    if saveFigs
        saveas(h1, strcat('venn-',tleg), 'png');
    end
end

%% Create a list of the hits overlapping between the different combinations of the time points
final_hits = {}; % will store the names of drug hits for each drug type
final_hits_OD = {}; % will store the OD of drug hits for each drug type
final_hits_metadata = {}; % will store the metadata of drug hits for each drug type
for iCond = 2:nConditions
    cur_tf_hits = tf_hits{1,iCond}; % extracts hits for condition i
%     cur_tf = sum([cur_tf_hits(:,inx_tp)],2)==2; % finds common hits at selected time points
    cur_tf = sum([cur_tf_hits(:,inx_tp)],2)>=1; % finds union of hits at selected time points
    [final_hits_OD{1,iCond-1}, orderHits] = sort(normOD(cur_tf,ODnormTime));
    % extracts OD of hits and sorts them 
    junk = drugArray(cur_tf); % intermediate variable holding unsorted drug names
    final_hits{1,iCond-1} = junk(orderHits); % stores drug names sorted by OD
    junk = Drugs(cur_tf,:);  % intermediate variable holding unsorted drug metadata
    final_hits_metadata{1,iCond-1} = junk(orderHits,:);  % stores drug metadata sorted by OD
end
% save filetosave final_hits final_hits_OD final_hits_metadata cutoffOD normOD ODnormTime tf_hits inx_conditions Drugs
%% Plots pie charts
h1 = figure; hold on;
pie([numel(inx_conditions{1,2})-numel(final_hits{1,1}) numel(final_hits{1,1})...
    numel(inx_conditions{1,3})-numel(final_hits{1,1}) numel(final_hits{1,1})])
axis off
title('Mega screen #3')
set(gcf,'position',[444,201,535,452])
legend({legendConditions{2} [legendConditions{2} ' hits'] legendConditions{3} [legendConditions{3} ' hits']},'location','southwest')

if saveFigs
    saveas(h1, 'pie_screened_and_hits', 'png');
end

%% Plots venn diagrams
h1 = figure; hold on;
for i=2:nConditions
    subplot(1,2,i-1)
    venn([numel(inx_conditions{1,i}),numel(final_hits{1,i-1})],numel(final_hits{1,i-1}))
    axis off
    text([-7 7],[1 1], {num2str(numel(inx_conditions{1,i})) num2str(numel(final_hits{1,i-1}))})
    title(legendConditions{i,1})
    set(gcf,'position',[440,465,927,333])
    legend({'screened' 'hits'})
end
if saveFigs
    saveas(h1, 'venn_screened_and_hits', 'png');
end
%% Plots barchart of hits
% barchart to have a quick look at ordered hits
for iCond = 1:nConditions-1
    h1 = figure;
    barh(final_hits_OD{1,iCond});
    set(gca,'YTick',[1:length(final_hits_OD{1,iCond})],'YTicklabel', final_hits{1,iCond});
    grid on;
    set(gcf, 'Position', [717,8,646,797]) % adjust the figure size for 2
% timepoints
%     set(gcf, 'PaperPosition', [0 0 20 20])% adjust the figure size for 12
% timepoints
    % title and file saving for antibiotic set
    title(strcat(legendConditions{iCond+1},' hits OD at-', num2str(ODnormTime-1), '-hours'))
    if saveFigs
        saveas(h1, strcat(legendConditions{iCond+1},'_hits_OD_',num2str(ODnormTime-1),'hr', '.png'));
    end
end
%% Saves xlsx file with all drugs' metadata sorted by normalized OD
antihost = Drugs(inx_conditions{1,ah},:);
[sortOD, inxOD] = sort(normOD(inx_conditions{1,ah},2));
antihost = antihost(inxOD,:);
antihost(:,end+1) = num2cell(sortOD);

writecell(antihost, xfiletosave)
