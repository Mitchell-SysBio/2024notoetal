%% Antibiotic screen in e. coli
% data analysis and plotting
% MN - November 08, 2021
% 384wp screen of del-ybcN (~WT e. coli) against 112 antibiotics (24 controls) 
% this script performs t-tests comparing the absorbance of drug treatments 
% to the absorbance when there's no drug
% Drugs with an fdr adjusted p-value smaller than 0.1 in any
% time point starting at 6 hours are considered a hit 
%% User definitions
saveFigs = false; % if true save figures created in analysis
% xfiletosave = strcat('20230316AntibioticHits','.xlsx');
nRep = 3;
nControl = 24;
myColormap = jet;
timeRes = 0.5; % in hours, measurements every 30 min
fdrThresh = 0.1;
%% load drug info
[~, drugNames] = xlsread('MCE_antibiotic_info.xlsx', '1-100_dilution', 'H2:H114'); 
[~, ~, drugID] = xlsread('MCE_antibiotic_info.xlsx', '1-100_dilution', 'A2:A114'); 
[~, ~, drugPos] = xlsread('MCE_antibiotic_info.xlsx', '1-100_dilution', 'E2:E114'); 
[~, ~, drugPlate] = xlsread('MCE_antibiotic_info.xlsx', '1-100_dilution', 'D2:D114'); 
drugArray = drugNames(1:end-1);
controlArray = drugNames(end);
%% load data
[rawPlateMap] = xlsread('MCE_antibiotic_info.xlsx','20211105exp','B2:Y17');
rawOD = xlsread('Mariana_OD600_20211106_104922.xlsx','D52:NW83');
[~,wellAddress] = xlsread('Mariana_OD600_20211106_104922.xlsx', 'D51:NW51'); %load well names
[~,wellIndex] = natsort(wellAddress); %get index to sort A1, A2,...,H12
OD600 = rawOD(:,wellIndex); %sort A1, A2,...,H12
inx384 = nan(16,24);
plateMap = [];
for i = 1:16
    inx384(i,:) = [(i-1)*24+1:i*24];
    plateMap = cat(2,plateMap,rawPlateMap(i,:));
end
OD600(:,isnan(plateMap)) = nan; % replace 
% this part of the code creates an array with the location of every drug
% and its replicates in the 384 well plate
nDrugs = numel(drugArray);
% create drug replicates index in 384 wp
inxDrugs = nan(nDrugs,nRep);
for i =1:nDrugs
    curDrugID = drugID{i};
    inxDrugs(i,:) = find(plateMap==curDrugID);
end
inxDMSO = find(plateMap==drugID{end});
nTime = size(OD600,1);
timeArray = [0:nTime-1].*timeRes;
OD600 = OD600-nanmean(OD600(1:3,:));

%% plate growth curves
% plots all replicates in one 384 well plate in the 112 format that
% summarizes the replicates
f = figure;
for iDrug = 1:nDrugs
    curInx = inxDrugs(iDrug,:);
    subplot(8,14,iDrug)
    plot(timeArray,OD600(:,curInx))
    xlabel('hours');
    ylabel('A600');
    title(['drug ' num2str(iDrug)])
    set(gca,'ylim',[min(min(OD600)) max(max(OD600))]);
    grid on; box on; axis square;
end
set(gcf,'position',[1 1 1440 804]);
if saveFigs
    saveas(f,'A600 overtime.png');
end
%%%%% Manually flag replicates that are off - order: blue; orange; yellow
% 2b 3b 12o 13b 16y 17y 18y 26y 28y 40b 54b 67o 70y 103y
inxMask = [2 1; 3 1; 12 2;13 1; 16 3; 17 3; 18 3; 26 3; 28 3; 40 1; 54 1;67 2; 70 3; 103 3];
for i = 1: size(inxMask,1)
    OD600(:,inxDrugs(inxMask(i,1),inxMask(i,2))) = nan;
end
f = figure;
for iDrug = 1:nDrugs
    curInx = inxDrugs(iDrug,:);
    subplot(8,14,iDrug); hold on;
    plot(timeArray,OD600(:,curInx))
    plot(timeArray,mean(OD600(:,inxDMSO),2),'-k','LineWidth',1)
    xlabel('hours');
    ylabel('A600');
    set(gca,'ylim',[min(min(OD600)) max(max(OD600))]);
    title(['drug ' num2str(iDrug)])
    grid on; box on; axis square;
end
sgtitle('growth curves - off replicates out');
set(gcf,'position',[1 1 1440 804]);
if saveFigs
    saveas(f,'A600 overtime - filtered.png');
end
%% Heatmap of 384 A600 values over time
% heatmap of raw absorbance reads in 384 well plates to visualize if there
% are any spatial patterns or biases
f = figure;
k=1;
for iTime = 1:2:nTime
    curmeanOD = OD600(iTime,:);
    plateOD = curmeanOD(inx384);
    subplot(4,4,k)
    heatmap(plateOD,'ColorLimits',[nanmedian(curmeanOD(:)) max(curmeanOD(:))],'colormap', myColormap)
    title(strcat(num2str(timeArray(iTime)), ' hours'))
    k = k+1;
end
sgtitle('heatmap A600 over time - 384wp format')
set(gcf,'position',[1 1 1440 804])
if saveFigs
    saveas(f, 'A600 384 wp over time.png');
end
%% Heatmap, t-test comparison and normalize OD by DMSO controls per column
% this section goes by drug and performs a t-test comparing the drug
% treated absorbance to the control in the same column. Specifically, it
% takes the value between the two replicates of the drug. Then, it estimates
% the positive false discovery rate for multiple hypothesis testing and
% makes a list of drugs with an adjusted p-value smaller than 0.25. This is
% done for all timepoints
% aditionally, it calculates the normalized OD to the same column control
p = nan(nTime,nDrugs); h = nan(nTime,nDrugs);fdr = nan(nTime,nDrugs); q = nan(nTime,nDrugs);
normOD = nan(nTime,nDrugs);
meanNormOD = nan(nTime,nDrugs);
meanOD = nan(nTime,nDrugs);

f = figure;
inxControl = inxDMSO;
k=1;
for iTime = 1:2:nTime
    controlOD = OD600(iTime,inxControl);
    curmeanOD = [];
    for iDrug = 1:nDrugs
        %             curInx = setdiff(inxDrugs(inxCurDrug,:),inxExclude);
        curInx = inxDrugs(iDrug,:);
        curOD = OD600(iTime,curInx);
        [h(iTime,iDrug), p(iTime,iDrug)] = ttest2(curOD,controlOD, 'tail','left','alpha', 0.05);
%         [h(iTime,iDrug), p(iTime,iDrug)] = ttest(curOD,controlOD(randperm(19,3)),'alpha',0.1, 'tail','left');
%         [h(iTime,iDrug), p(iTime,iDrug)] = ttest(curOD,controlOD(randperm(19,3)),'alpha',0.1, 'tail','left');
        normOD(iTime,curInx) = curOD./nanmean(controlOD); %normalize OD by mean OD of all controls
        meanNormOD(iTime,iDrug) = nanmean(normOD(iTime,curInx));
        curmeanOD = [curmeanOD,nanmean(curOD)];
    end
%     [fdr(iTime,:),q(iTime,:)] = mafdr(p(iTime,:),'bhfdr',1);
    [fdr(iTime,:)] = mafdr(p(iTime,:),'bhfdr',1);
    plateMeanOD = reshape(curmeanOD,[14,8])';
    plateMeanOD = plateMeanOD([8:-1:1],:);
    subplot(4,4,k); hold on;
    imagesc(plateMeanOD,[min(plateMeanOD(:)) max(plateMeanOD(:))])%,'colormap', myColormap)
    colorbar;
    axis tight; colormap(myColormap)
    l = 8;
    for iRow = 1:8
        for iCol = 1:14
            %                 curh = h(iTime,iCol+(iRow-1)*12+(iPlate-1)*96);
            %                 if curh, plot(iCol,iRow,'.k'); end
            curFDR = fdr(iTime,iCol+(iRow-1)*14);
            if curFDR < fdrThresh, plot(iCol,l,'.k'); end
        end
        l=l-1;
    end
    title(strcat(num2str(timeArray(iTime)),' hours'))
    k = k+1;
    meanOD(iTime,:) = curmeanOD;
end
sgtitle('A600 over time')
set(gcf,'position',[1 1 1440 804])
if saveFigs
    saveas(f,'heatmap A600 over time with hits.png');
end
f = figure;
k=1;
for iTime = 1:2:nTime
    plateMeanNormOD = reshape(meanNormOD(iTime,:),[14,8])';
    plateMeanNormOD = plateMeanNormOD([8:-1:1],:);
    subplot(4,4,k); hold on;
    imagesc(plateMeanNormOD,[min(plateMeanNormOD(:)) max(plateMeanNormOD(:))])%,'colormap', myColormap)
    axis tight; colormap(myColormap); colorbar;
    l = 8;
    for iRow = 1:8
        for iCol = 1:14
            %                 curh = h(iTime,iCol+(iRow-1)*12+(iPlate-1)*96);
            %                 if curh, plot(iCol,iRow,'.k'); end
            curFDR = fdr(iTime,iCol+(iRow-1)*14);
            if curFDR < fdrThresh, plot(iCol,l,'.k'); end
        end
        l = l - 1;
    end
    title(strcat(num2str(timeArray(iTime)),' hours'))
    k = k+1;
end
sgtitle('Norm A600 over time')
set(gcf,'position',[1 1 1440 804])
if saveFigs
    saveas(f,'heatmap A600 over time normalized to control with hits.png');
end
%% Make table with hits through time
% this section finds drugs that are hits every hour starting at 5 hours 
% and extracts the normalized OD and metadata for each hit
tf = fdr<fdrThresh;
t = [11:2:nTime];
inxHits = tf(t,:); % hit in at least 2 time pointx
nHits = [(t'-1)/2 sum(inxHits,2)];
f=figure; hold on;
bar((t'-1)/2,nHits(:,2));
yline(112,'--r','lineWidth',2);
text(12, 110,'screened antibiotics')
xlabel('hours')
ylabel(['number of hits - FDR' num2str(fdrThresh)])
grid on; box on;
if saveFigs
    saveas(f,'number of hits over time.png');
end
%% make histograms of OD and normalized OD over time
edges = [0.15:0.1:1.25];
f = figure;
k=1;
for iTime = 1:2:nTime
    subplot(4,4,k)
    histogram(meanNormOD(iTime,:),edges,'normalization','probability')
    box on; grid on;
    xlabel('normalized A600');
    ylabel('proportion');
    title([num2str(timeArray(iTime)),' hours - ',num2str(sum(tf(iTime,:))), ' hits'])
    k = k + 1;
end
set(gcf,'position',[382 105 1011 597]);
sgtitle('histogram normalized A600 over time');
if saveFigs
    saveas(f,'histogram normalized A600 over time.png');
end

edges = [0.09:0.04:1.1];
f = figure;
k=1;
for iTime = 1:2:nTime
    subplot(4,4,k); hold on;
    histogram(meanOD(iTime,:),edges,'normalization','probability')
    histogram([nanmean(OD600(iTime,inxDMSO)) 2],edges,'normalization','probability','facecolor','r')
    box on; grid on;
    xlabel('A600');
    ylabel('proportion');
    title([num2str(timeArray(iTime)),' hours - ',num2str(sum(tf(iTime,:))), ' hits'])
    k = k + 1;
end
set(gcf,'position',[382 105 1011 597]);
sgtitle('histogram A600 over time');
if saveFigs
    saveas(f,'histogram A600 over time.png');
end

%% save data for supp table 1 2024-01-12
t = [19:2:nTime];
hitsnormOD = nan(length(t),nDrugs);
hitsFDR =nan(length(t),nDrugs);
hitsNames = {};
for i=1:length(t)
    hitsnormOD(i,:) = meanNormOD(t(i),:);
    hitsFDR(i,:) = fdr(t(i),:);
    hitsNames(:,i) = drugNames(1:end-1);
    T = table(hitsNames(:,i),...
    hitsnormOD(i,:)',hitsFDR(i,:)',...
    'VariableNames',{'drug_names', 'normA600','FDR'});
    writetable(T,[date 'allNormOD.xlsx'],'sheet',['time_', num2str((t(i)-1)/2), '_hours'])  
end