%% Subset of FDA library drug screen in e. coli
% data analysis and plotting
% June 06, 2021
% 2X384wp screen of del-ybcN (~WT e. coli) against 176 drugs (16 controls)
% this script joins data from different plates and runs t-tests comparing
% the absorbance of drug treatments to the fit of the absorbance of the
% control wells taking into account the location in the plate, specifically
% the column. Drugs with an fdr adjusted p-value smaller than 0.25 in any
% time point starting at 6 hours are considered a hit and will be selected
% to be included in the ASKA screen
%% User definitions
saveFigs = false; % if true save figures created in analysis
xfiletosave = strcat(date ,'_AntihostHits_FDR','.xlsx');
%% load data
nControls = 8; %number of controls per plate, not considering duplicates
nPlates = 2; % number of 384 well plates
myColormap = jet;

% this part of the screen loads the drug nameas and sorts A1, A2,...,B1,...,H12
% using the well names provided in the file
[~, drugNames] = xlsread('20210603sortedAntihostDrugs.xlsx', 'Sheet3', 'B1:I96'); %load drug names for plate 1
[~, wellAddress] = xlsread('20210603sortedAntihostDrugs.xlsx', 'Sheet3', 'F1:F96'); %load drug wells
[drugWells,wellIndex] = natsort(wellAddress); %get index to sort A1, A2,...,H12
drugNames1 = drugNames(wellIndex,1); %sort A1, A2,...,H12
junk = [2:3 6:8];
drugMeta1 = drugNames(wellIndex,junk);
[~, drugNames] = xlsread('20210603sortedAntihostDrugs.xlsx', 'Sheet3', 'B97:I192'); %load drug names for plate 1
[~, wellAddress] = xlsread('20210603sortedAntihostDrugs.xlsx', 'Sheet3', 'F97:F192');  %load drug wells
[~,wellIndex] = natsort(wellAddress); %get index to sort A1, A2,...,H12
drugNames2 = drugNames(wellIndex,1);  %sort A1, A2,...,H12
drugMeta2 = drugNames(wellIndex,junk);

drugArray = {drugNames1{:}, drugNames2{:}}; % make drug array of all plates
drugMeta = {drugMeta1{:,1}, drugMeta2{:,1};drugMeta1{:,2}, drugMeta2{:,2};...
    drugMeta1{:,3}, drugMeta2{:,3};drugMeta1{:,4}, drugMeta2{:,4};drugMeta1{:,5}, drugMeta2{:,5}}';

% this part of the code creates an array with the location of every drug
% and its replicates in the 384 well plates
nDrugs = numel(drugArray);
% create drug replicates index in 384 wp
inxDrugs = [1:2];
for i =1:nDrugs/2-1
    inxDrugs = [inxDrugs; inxDrugs(end,:)+2];
end
inxDrugs = [inxDrugs, inxDrugs+(96*2)];
inxDrugs = [inxDrugs; inxDrugs+384];

% this part of the screen loads the data and sorts A1, A2,...,B1,...,H12
% using the well names provided in the Tecan output file
rawOD = xlsread('OD_20210605_plate1.xlsx', 'D51:NW63'); %load data from plate1
[~,wellAddress] = xlsread('OD_ecoli_20210605_plate1.xlsx', 'D50:NW50'); %load well names
[~,wellIndex] = natsort(wellAddress);
plate1 = rawOD(:,wellIndex); %sort A1, A2,...,H12

rawOD = xlsread('OD_20210605_plate2.xlsx', 'D52:NW64'); %load data from plate2
[~,wellAddress] = xlsread('OD_20210605_plate2.xlsx', 'D51:NW51');
[~,wellIndex] = natsort(wellAddress);
plate2 = rawOD(:,wellIndex); %sort A1, A2,...,H12

OD600 = [plate1(2:end,:)-plate1(1,:),plate2(2:end,:)-plate2(1,:)]; %create
%1 matrix with all OD values and delete values from timepoint zero after background subtraction

clear rawOD plate1 plate2 wellAddress wellIndex drugMeta1 drugMeta2 junk
nTime = size(OD600,1);
timeArray = [1:nTime];

%% plate growth curves
% plots all replicates in one 384 well plate in the 96 format that
% summarizes the replicates
for iPlate = 1:nPlates
    f = figure;
    for iDrug = 1:nDrugs/nPlates
        curInx = inxDrugs((iPlate-1)*nDrugs/2+iDrug,:);
        subplot(8,12,iDrug)
        plot(timeArray,OD600(:,curInx))
        set(gca,'ylim',[min(min(OD600)) max(max(OD600))])
        title(drugArray{(iPlate-1)*nDrugs/2+iDrug})
        grid on; box on; axis square;
    end
    sgtitle(strcat('A600 - Plate',num2str(iPlate)))
    set(gcf,'position',[12,78,1392,726])
    if saveFigs
        saveas(f, strcat('A600 - Plate',num2str(iPlate)), 'png');
    end
end
%% Heatmap of 384 A600 values over time
% heatmap of raw absorbance reads in 384 well plates to visualize if there
% are any spatial patterns or biases
for iPlate = 1:nPlates
    f = figure;
    if(iPlate == 1), curInx = 1:384;
    else, curInx = 385:768; end

    curPlate = OD600(:,curInx);
    for iTime = 1:nTime
        curmeanOD = curPlate(iTime,:);
        plateOD = reshape(curmeanOD,[24,16])';
        subplot(3,4,iTime)
        heatmap(plateOD,'ColorLimits',[median(curmeanOD(:)) max(curmeanOD(:))],'colormap', myColormap)
        title(strcat('hour',num2str(timeArray(iTime))))
    end
    sgtitle(strcat('A600 over time - Plate',num2str(iPlate)))
    set(gcf,'position',[90,48,1256,549])
    if saveFigs
        saveas(f, ['A600 384 wp over time - Plate',num2str(iPlate),' (median to max)','.png']);
    end
end
%% Build controls
% takes control and fits a model to tha absorbance as a function of
% position in plate. This is done in order to normalize absorbance to the
% controls in the same column of the plate, to account for small position
% biases seen on the 384 heatmaps
inxDMSO = inxDrugs(find(strcmp(drugArray,'DMSO')),:);
x =  mod(inxDMSO,24);
F = {};
for iPlate = 1:nPlates
    inxControl = intersect(inxDMSO(:), [1:384]+384*(iPlate-1));
    curX = x(1:8,1:4);
    curY = inxDMSO(1:8,1:4)+(iPlate-1)*384;
    for iTime = 1:nTime
        F{iTime,iPlate} = fit(curX(:),OD600(iTime,curY(:))','c1*(exp(-x*c2))+c3','StartPoint',[range(OD600(iTime,curY(:))) 0.05 min(OD600(iTime,curY(:)))],'MaxIter',1000);
    end
end
%% plot controls and fit function
% visualize controls and fit
colormap(jet)
f= figure;
for iPlate = 1:nPlates
    for iTime =1:nTime
        %     plot(sort(p(iTime,:)))
        subplot(3,4,iTime); hold on;
        curX = x(1:8,1:4);
        curY = inxDMSO(1:8,1:4)+(iPlate-1)*384;
        plot(curX(:),OD600(iTime,curY(:)),'.')
        plot([1:0.1:24],F{iTime,iPlate}([1:0.1:24]))
        set(gca,'xlim',[1 24])
        legend({'plate1','fit1','plate2','fit2'})
        xlabel('column position')
        ylabel('DMSO A600')
    end
    set(gcf, 'position', [1,194,999,604]);
    if saveFigs
        saveas(f, ['DMSO controls over time ','.png']);
    end
end
%% Heatmap, t-test comparison and normalize OD by DMSO controls per column
% this section goes by drug and performs a t-test comparing the drug
% treated absorbance to the control in the same column. Specifically, it
% takes the value between the two replicates of the drug. Then, it estimates
% the positive false discovery rate for multiple hypothesis testing and
% makes a list of drugs with an adjusted p-value smaller than 0.25. This is
% done for all timepoints
% aditionally, it calculates the normalized OD to the same column control
inxExclude = [1:24:384];
p = nan(nTime,nDrugs); h = nan(nTime,nDrugs);fdr = nan(nTime,nDrugs); q = nan(nTime,nDrugs);
normOD = nan(nTime,nDrugs);
meanNormOD = nan(nTime,nDrugs/nPlates);
for iPlate = 1:nPlates
    f = figure;
    inxControl = intersect(inxDMSO(:), [1:384]+384*(iPlate-1));
    for iTime = 1:nTime
        controlOD = OD600(iTime,inxControl);
        curmeanOD = [];
        for iDrug = 1:nDrugs/nPlates
            inxCurDrug = (iPlate-1)*nDrugs/2+iDrug; % goes from 1 to 192
            %             curInx = setdiff(inxDrugs(inxCurDrug,:),inxExclude);
            curInx = inxDrugs(inxCurDrug,:);
            curOD = OD600(iTime,curInx);
            [h(iTime,inxCurDrug), p(iTime,inxCurDrug)] = ttest(curOD,F{iTime,iPlate}(mod(min(curInx)+0.5,24)),'alpha',0.25, 'tail','left');
            normOD(iTime,curInx) = curOD./F{iTime,iPlate}(mod(min(curInx)+0.5,24)); %normalize OD by fitted value from the midpoint of the two consecutive replicates
            meanNormOD(iTime,inxCurDrug) = mean(normOD(iTime,curInx));
            curmeanOD = [curmeanOD,mean(curOD)];
        end
        %         [fdr(iTime,:),q(iTime,:)] = mafdr(p(iTime,:));
        fdr(iTime,:) = mafdr(p(iTime,:),'bhfdr',1);
        plateMeanOD = reshape(curmeanOD,[12,8])';
        subplot(3,4,iTime); hold on;
        imagesc(plateMeanOD,[min(plateMeanOD(:)) max(plateMeanOD(:))])%,'colormap', myColormap)
        axis tight; colormap(myColormap)
        for iRow = 1:8
            for iCol = 1:12
                %                 curh = h(iTime,iCol+(iRow-1)*12+(iPlate-1)*96);
                %                 if curh, plot(iCol,iRow,'.k'); end
                curFDR = fdr(iTime,iCol+(iRow-1)*12+(iPlate-1)*96);
                if curFDR<0.25, plot(iCol,iRow,'.k'); end
            end
        end
        title(strcat('hour',num2str(timeArray(iTime))))
    end
    sgtitle(strcat('A600 over time - Plate',num2str(iPlate)))
    set(gcf,'position',[90,48,1256,549])
    if saveFigs
        saveas(f, strcat('A600 over time - Plate',num2str(iPlate)), 'png');
    end
end

%% Make table with hits for timepoint T
% this section finds drugs that are hits in any of the timepoints after 6
% hours and extracts the normalized OD and metadata for each hit
tf = fdr<0.25;
t = [6:1:nTime];
T=9;
inxHits = find(sum(tf(t,:))>1); % hit in at least 2 time pointx
curHits = drugArray(inxHits)';
hitsMeanNormOD = meanNormOD(t,inxHits)';
hitsFDR = fdr(t,inxHits)';
hitsMeta = drugMeta(inxHits,:);
junk = [drugWells,drugWells];
hitsWell = junk(inxHits);
[~, inxSortFDR] = sort(fdr(T,inxHits)); %sort by FDR at timepoint T
%%
hitTable = table(curHits(inxSortFDR),hitsMeta(inxSortFDR,:),hitsWell(inxSortFDR)',hitsFDR(inxSortFDR,:),hitsMeanNormOD(inxSortFDR,:),'VariableNames',{'Antihost drugs','Metadata','position in 96',['fdr t', num2str(t(1)),' to t',num2str(t(end)),' sorted at t', num2str(T)],['mean normalized OD t', num2str(t(1)),' to t',num2str(t(end)),' sorted at t', num2str(T)]});
% hitTable = table(curHits(inxSortFDR),hitsFDR(inxSortFDR,1),hitsFDR(inxSortFDR,2),hitsFDR(inxSortFDR,3),hitsFDR(inxSortFDR,4),hitsFDR(inxSortFDR,5),hitsMeanNormOD(inxSortFDR,1),hitsMeanNormOD(inxSortFDR,2),hitsMeanNormOD(inxSortFDR,3),hitsMeanNormOD(inxSortFDR,4),hitsMeanNormOD(inxSortFDR,5),'VariableNames',{'Drug','fdr t=8','fdr t=9','fdr t=10','fdr t=11','fdr t=12','normalized OD t=8','normalized OD t=9','normalized OD t=10','normalized OD t=11','normalized OD t=12'});
writetable(hitTable,xfiletosave)
%%
save miniscreen20210604_0617.mat drugArray meanNormOD inxHits curHits hitsMeanNormOD hitsFDR fdr inxSortFDR
%% pie chart with drugs and hits
f = figure;
pie([nDrugs-(nControls*nPlates)-numel(curHits) numel(curHits)])
legend({'antihost drugs','antihost drug hits'},'Location','northwest')
title('Top 172 antihost drugs')
if saveFigs
    saveas(f, 'pieAntihostHits','png');
end
%% plot the histogram of normOD by time point
h = figure; hold on;
edges = [-0.01:0.02:1.11];
inx_tp = [7:nTime]; % time points to inspect
n = length(inx_tp); % useful variable for later


for i=1:n
    cur_inx_tp = inx_tp(i); % the the index of the time point for the current iteration
    subplot(1,n,i); hold on;
    histogram(meanNormOD(cur_inx_tp,:),edges,'normalization','probability');
    xline(max(meanNormOD(cur_inx_tp,fdr(cur_inx_tp,:)<0.25)),'k'); % plot a line to show the cutoff
    n_drugs = nDrugs-(nControls*nPlates); % the number of wells in the distribution
    cur_n_hits = sum(fdr(cur_inx_tp,:)<0.25); % number of toxic drug by the cutoff
    % we build a string with the values we just calculated
    title({['hour ' num2str(cur_inx_tp) ' (', num2str(cur_n_hits), '/', num2str(n_drugs), ')']});
    box on; axis square;
end
sgtitle('normalized OD')
set(h,'Position',[5,433,1436,188]); % adjust the figure size for 2
if saveFigs
    saveas(h,'normalized_OD_by_drug_type_decreased_growth_hits', 'png');
end

%% plot the histogram of fdr by time point
h = figure; hold on;
edges = [-0.01:0.02:0.73];
inx_tp = [7:nTime]; % time points to inspect
n = length(inx_tp); % useful variable for later


for i=1:n
    cur_inx_tp = inx_tp(i); % the the index of the time point for the current iteration
    subplot(1,n,i); hold on;
    histogram(fdr(cur_inx_tp,:),edges,'normalization','probability');
    xline((0.25),'k'); % plot a line to show the cutoff
    n_drugs = nDrugs-(nControls*nPlates); % the number of wells in the distribution
    cur_n_hits = sum(fdr(cur_inx_tp,:)<0.25); % number of toxic drug by the cutoff
    % we build a string with the values we just calculated
    title({['hour ' num2str(cur_inx_tp) ' (', num2str(cur_n_hits), '/', num2str(n_drugs), ')']});
    box on; axis square;
end
sgtitle('fdr adjusted p-values')

set(h,'Position',[5,433,1436,188]); % adjust the figure size for 2
if saveFigs
    saveas(h,'fdr_adjusted_pvalues', 'png');
end

