%% Efflux systems analysis
% this script generates analysis and plots for figure 4
% 11/26/2022
%% Data preparation
cutoffpadj = 0.25; % threshold to set  hits based on p-val
% Load annotation data
load strainResultsAggregated.mat; %
raw = readcell('metadata/drugMetadata.xlsx');
strainMetadata = readcell('metadata/strainMetadata.csv');
% Load screen data
% The raw drug metadata file is not in the same order as the screens in
% strainResults and includes ignored screens so we need to sort it first
% and then filter
screenID = strainResults.screenID'; % gets screen IDs
% screenID(45) = []; disp('Message: Screen #45 was manually removed')
nScreens = length(screenID); % useful variable

screenInx2 = zeros(length(screenID),1);
for i=1:nScreens
    curID = screenID{i};
    if ismember(i,166:167)
        screenInx2(i) = find(contains(raw(:,strcmp(raw(1,:),'fileRep2')),curID));
    else
        if isempty(find(contains(raw(:,strcmp(raw(1,:),'fileRep1')),curID)))
            continue
        else
            screenInx2(i) = find(contains(raw(:,strcmp(raw(1,:),'fileRep1')),curID));
        end
    end
end

strainNames = strainResults.strainNames;
drugMetadata = raw(screenInx2,:); % get medatadata on drugs sorted & filtered
drugMetadataHeader = raw(1,:);

% Get LFC and LFC masked by padj
LFC = strainResults.LFC; % LFC of filtered screens
padj = strainResults.padj; % padj of filtered screens
TFpadj = padj < cutoffpadj; % TF of screens with padj over threshold
LFCmasked = strainResults.LFCmasked;
padjMasked = nan(size(padj)); % replace padh of non-hits with
padjMasked(TFpadj) = padj(TFpadj);
hitsMat = ~isnan(LFCmasked);
% standardize screens
% LFC(:,45)=[]; disp('Message: Screen #45 was removed from LFC')
% TFpadj(:,45)=[]; disp('Message: Screen #45 was removed from TFpadj')
stdLFC = nanstd(LFC,0,1);
normLFC = LFC./stdLFC;
normLFCmasked = normLFC;
normLFCmasked(~TFpadj) = 0;
data = normLFC';
imputedData = knnimpute(data,5);
normLFC2 = imputedData';
data = normLFCmasked';
imputedData = knnimpute(data,5);
normLFCmasked2 = imputedData';
% Parse annotation data
mechTerm = 'mechanism'; % 'mechanism' for targeted mechanism, 'class' for chemical class
inxMech = strcmp(drugMetadataHeader,mechTerm);

% get TF where 1 means the screen was with an antibiotic
tfUnknown = logical(strcmp(drugMetadata(:,inxMech),'unknown')); % drugs with unknown class/mechanism
tfHt = logical(cell2mat(drugMetadata(:,strcmp(drugMetadataHeader,'antihost')))); % only host-targeted drugs
tfAbAll = logical(cell2mat(drugMetadata(:,strcmp(drugMetadataHeader,'antibiotic')))); % include unclassified antibiotics
tfAb = tfAbAll & ~tfUnknown;

drugNames = drugMetadata(:,strcmp(drugMetadataHeader,'Name')); % array to label with drug names
drugMechanismRaw = drugMetadata(:,strcmp(drugMetadataHeader,mechTerm));
uniqMechanism = unique(drugMechanismRaw(tfAbAll)); % unique mechanism terms
uniqMechanism{end+1} = 'anti-host';
ATCraw = drugMetadata(:,strcmp(drugMetadataHeader,'ATC1'));
% Parse mechanisms and prepare mechanism arrays
mechanismPrimary = cell(nScreens,1); mechanismPrimaryInx = nan(nScreens,1);
mechanismSecondary = cell(nScreens,1); mechanismSecondaryInx = nan(nScreens,1);
mechanismPrimaryCounts = zeros(size(uniqMechanism));

for i = 1:nScreens
    curMech = drugMechanismRaw{i};
    if(tfHt(i)) % antihost drug
        mechanismPrimary{i} = 'anti-host';
        mechanismSecondary{i} = drugMechanismRaw{i}; mechanismSecondaryInx(i) = 3;
    else % antibiotic drug
        if(tfAb(i)) % antibiotic with known mechanism
            mechanismPrimary{i} = drugMechanismRaw{i};
            mechanismSecondary{i} = drugMechanismRaw{i}; mechanismSecondaryInx(i) = 1;
        else % antibiotic with unknown mechanism
            mechanismPrimary{i} = drugMechanismRaw{i};
            mechanismSecondary{i} = 'ab-unknown'; mechanismSecondaryInx(i) = 2;

        end
    end
    curPrimaryInx = find(strcmp(uniqMechanism,mechanismPrimary{i}));
    mechanismPrimaryInx(i) = curPrimaryInx;
    mechanismPrimaryCounts(curPrimaryInx) = mechanismPrimaryCounts(curPrimaryInx)+1;
end

% house cleaning
clear raw i curID inxMech inxName;
clear drugMechanismRaw uniqAbMechanism curPrimaryInx;
% Parse HT mechanisms and prepare mechanism arrays
drugMechanismRaw = drugMetadata(:,strcmp(drugMetadataHeader,mechTerm));
uniqMechanismHt = unique(drugMechanismRaw(tfHt)); % unique mechanism terms
uniqMechanismHt{end+1} = 'ab';
uniqMechanismHt{end+1} = 'ab-unknown';

mechanismHt = cell(nScreens,1); mechanismHtInx = nan(nScreens,1);
mechanismHtCounts = zeros(size(uniqMechanismHt));

for i = 1:nScreens
    curMech = drugMechanismRaw{i};
    if(tfHt(i)) % antihost drug
        mechanismHt{i} = drugMechanismRaw{i};
    else % antibiotic drug
        if(tfAb(i)) % antibiotic with known mechanism
            mechanismHt{i} = 'ab';
        else % antibiotic with unknown mechanism
            mechanismHt{i} = 'ab-unknown';
        end
    end
    curPrimaryInx = find(strcmp(uniqMechanismHt,mechanismHt{i}));
    mechanismHtInx(i) = curPrimaryInx;
    mechanismHtCounts(curPrimaryInx) = mechanismHtCounts(curPrimaryInx)+1;
end
% clear drugMechanismRaw curPrimaryInx;
% Build colormap
mechColor = [186,85,211; 125,25,207; 220,20,60; 255,140,0; 255,215,0; 0,100,100; 0,128,0; 0,0,180; 240,128,128; 218,112,214; 100,100,100];
% mechColor = mechColor./max(mechColor,[],2);
mechColor = mechColor./255;
mechColor(end,:) =  [0.7 0.7 0.7];
mechColor(end+1,:) =  [0.3 0.3 0.3];
drugColorLabel = zeros(length(mechanismPrimaryInx),3);
for i=1:size(drugColorLabel,1)
    drugColorLabel(i,:) = mechColor(mechanismPrimaryInx(i),:);
end

fh=figure('color','white');
p=pie(mechanismPrimaryCounts,[zeros(1,11) 1],num2str(mechanismPrimaryCounts));
legend(uniqMechanism,'location','westoutside');
colormap(mechColor);
% Infer useful indexes for drug groups and gene groups
growthRate= cell2mat(strainMetadata(ismember(strainMetadata(:,1),strainNames),17));
slowGrowerNames = strainNames(growthRate<0.7);
inxSlowStrains = find(ismember(strainNames,slowGrowerNames));

load ./metadata/prophageStrains.mat
prophages(strcmp(prophages(:,1),'yagN'),:)=[];

inxProphageStrains = find(ismember(strainNames,prophages(:,1)));

load ./metadata/pumpStrains.mat
inxPumpStrains = find(ismember(strainNames,pumps(:,1)));

inxMembraneStrains = unique([find(startsWith(strainNames,'omp'));find(startsWith(strainNames,'tol'));...
    find(startsWith(strainNames,'ton'));find(startsWith(strainNames,'pal'));find(startsWith(strainNames,'smpA'))...
    ;find(startsWith(strainNames,'yfgL'));find(startsWith(strainNames,'asmA'));find(startsWith(strainNames,'surA'));...
    find(startsWith(strainNames,'exbB'));find(startsWith(strainNames,'exbD'));find(startsWith(strainNames,'marR'));...
    find(startsWith(strainNames,'rob'));find(startsWith(strainNames,'hlpA')); find(startsWith(strainNames,'yrb'))]);

inxAllDrugs = [1:178]';
inxAbDrugs = find(tfAb);
inxHtDrugs = find(tfHt);
inxAllCharDrugs = [inxAbDrugs; inxHtDrugs];

inxAllStrains = 1:size(normLFC,1);
inxRemoveStrains = [inxSlowStrains;inxProphageStrains;inxPumpStrains;inxMembraneStrains];
inxAllUsefullStrains = setdiff(inxAllStrains,inxRemoveStrains);

TFLFC = abs(normLFC) > log2(2.5);
tfHits = TFLFC & TFpadj;
tfHitsPosNeg = double(tfHits);
tfHitsPosNeg(normLFC < -log2(2.5)) = 0.000005;
tfHitsPosNeg(tfHitsPosNeg==0) = 0.5;
% Ab based analysis
nStrainHits = sum(tfHits(:,inxAbDrugs),2);
inxAbInformativeStrains = find(nStrainHits>3 & nStrainHits <30);
inxStrains = setdiff(inxAbInformativeStrains,inxRemoveStrains);

%% Pump enrichment analysis
pumpsRaw = readcell(['pumpsManuallyMN_final.xlsx']);
pumpNames = unique(pumpsRaw(2:end,strcmp(pumpsRaw(1,:),'pump name')));
nPumps = numel(pumpNames);
tfNparts = zeros(nPumps,1);
singleStrains={}; 
for iPump=1:nPumps
    curName = pumpNames{iPump};
    curParts = pumpsRaw(strcmp(pumpsRaw(:,strcmp(pumpsRaw(1,:),'pump name')),curName),strcmp(pumpsRaw(1,:),'strain'));
    inxStrains = find(ismember(strainNames,curParts));
    pumpFam(iPump) = unique(pumpsRaw(strcmp(pumpsRaw(:,strcmp(pumpsRaw(1,:),'pump name')),curName),strcmp(pumpsRaw(1,:),'family')));
%      calculate pval (drugs vs. specific pump)
    for iDrug=1:length(inxAllDrugs)
        pval(iPump,iDrug) = calcEmpiricalPVal((normLFC2(:,inxAllDrugs(iDrug))),(normLFC2(inxStrains,inxAllDrugs(iDrug))));
    end
    if (numel(curParts)==1)
        tfNparts(iPump,1) = 1;
        singleStrains = cat(1,singleStrains,curParts);
        pval(iPump,:) = tfHitsPosNeg(inxStrains,:);
    end
%         fh = plotSingleHCluster(normLFC2(inxStrains,inxAllDrugs),strainNames(inxStrains),drugNames(inxAllDrugs),drugColorLabel(inxAllDrugs,:),curName,true);
%     saveas(fh,[curName '_hcluster.svg']); close;
%     fh = plotSingleHCluster(normLFCmasked2(inxStrains,inxAllDrugs),strainNames(inxStrains),drugNames(inxAllDrugs),drugColorLabel(inxAllDrugs,:),curName,true);
%     saveas(fh,[curName 'masked_hcluster.svg']); close;
end
% get pump hits
pumpImpact = zeros(size(pval));
pumpImpact(pval<0.01) = -1;
pumpImpact(pval>0.99) = 1;
% pumpImpactMasked = zeros(size(pvalMasked));
% pumpImpactMasked(pvalMasked<0.01) = -1;
% pumpImpactMasked(pvalMasked>0.99) = 1;

% figure; subplot(1,2,1);
bar(sum(pumpImpact~=0,2))
% subplot(1,2,2);
% bar(sum(pumpImpactMasked~=0,2))
% 
% inxStrains = ismember(strainNames,singleStrains);
% figure; plotSingleHCluster(normLFCmasked2(inxStrains,inxAllDrugs),strainNames(inxStrains),drugNames(inxAllDrugs),drugColorLabel(inxAllDrugs,:),'single strain pumps',true);
%% Pies - sensitivity increased by at least one pump
abSens = sum(sum(pumpImpact(:,tfAb)<0)>0); abTot = numel(inxAbDrugs);
htSens = sum(sum(pumpImpact(:,tfHt)<0)>0); htTot = numel(inxHtDrugs);
figure('Color','white');
subplot(1,2,1); pie([abTot-abSens abSens]);
subplot(1,2,2); pie([htTot-htSens htSens]);
%% pumps by family - sensitivity
nMech = numel(uniqMechanism);
famArray = unique(pumpFam);

plotInx = 1:numel(famArray);
f=figure('color', 'white'); t=tiledlayout(2,4);
sizePies = [];    X = zeros(numel(famArray),nMech);
for j =1:numel(plotInx)
    curFam = famArray{j};
    i=find(strcmp(pumpFam,curFam));
    %     curPump = pumpNames{i};
    curImpact = pumpImpact(i,:);
    for iMech = 1:nMech
        if length(i)>1
            X(j,iMech) = sum(sum(curImpact(:,mechanismPrimaryInx==iMech)<0)>0);
        else
            X(j,iMech) = sum((curImpact(:,mechanismPrimaryInx==iMech)<0)>0);
        end
    end
    nexttile; p=pie(X(j,:),string(X(j,:)));k=1;
    for iPie = 1:2:length(p)
        p(iPie).FaceColor = mechColor(k,:);
        k=k+1;
    end
    sizePies = [sizePies, sum(X(j,:))];
    title(curFam);
end
legend(uniqMechanism,'Location','northeast');
set(gcf,"Position",[1 83 1035 783]); title(t,'sensitivity')
% saveas(f,'piesPumpFamByMech.svg')

sizeFactor = sizePies*400;
sizeDiam = sqrt(sizeFactor/pi)*2;
sizeTable = table(famArray', sizePies', sizeFactor',sizeDiam');
%% contingency table and enrichment analysis
fisherP = zeros(numel(famArray),nMech);
fisherH = zeros(numel(famArray),nMech);
fdrP = zeros(numel(famArray),nMech);
for j =1:numel(famArray)
    curX = X(j,:);
    for iMech = 1:nMech
        curSens = curX(iMech);
        curTot = mechanismPrimaryCounts(iMech);
        ct = [curSens curTot-curSens; sum(curX)-curSens nScreens-curTot-(sum(curX)-curSens)];
        [fisherH(j,iMech), fisherP(j,iMech)] = fishertest(ct);
    end
    fdrP(j,:) = mafdr(fisherP(j,:),'bhfdr',1);

end
fdrH = fdrP<=0.055;
Th = table(famArray',fisherH);
Th = splitvars(Th,"fisherH",'NewVariableNames',uniqMechanism);
Tp = table(famArray',fisherP);
Tp = splitvars(Tp,"fisherP",'NewVariableNames',uniqMechanism);
Tf = table(famArray',fdrP);
Tf = splitvars(Tf,"fdrP",'NewVariableNames',uniqMechanism);
Tfh = table(famArray',fdrH);
Tfh = splitvars(Tfh,"fdrH",'NewVariableNames',uniqMechanism);
% writetable(Th,'pumpFisherComparison.xlsx','sheet','fisherH');
% writetable(Tp,'pumpFisherComparison.xlsx','sheet','fisherP');
% writetable(Tfh,'pumpFisherComparison.xlsx','sheet','fdrH');
% writetable(Tf,'pumpFisherComparison.xlsx','sheet','fdr');
%% pump sensitivity plot
[~,pumpImpactInx] = sort(sum(pumpImpact<0,2),'descend');
sortedPumpImpact = pumpImpact(pumpImpactInx,:);
sortedPumpNames = pumpNames(pumpImpactInx);
famSorted = pumpFam(pumpImpactInx);
[~,pumpImpactInx] = sort(famSorted);
sortedPumpImpact = sortedPumpImpact(pumpImpactInx,:);
sortedPumpImpact(sortedPumpImpact==1) = 0;
sortedPumpNames = sortedPumpNames(pumpImpactInx);
pumpImpactInx = [24:26 1:6 7 29:30 8:23 31 33 32];
sortedPumpImpact = sortedPumpImpact(pumpImpactInx,:);
sortedPumpNames = sortedPumpNames(pumpImpactInx);
nPumpByDrug = sum(pumpImpact<0);
% delPumps = find(sum(sortedPumpImpact~=0,2)==0);
% sortedPumpImpact(delPumps,:)=[];
% sortedPumpNames(:,delPumps)=[];
nPumps = numel(sortedPumpNames);
f1=figure('color', 'white'); hold on; lastN=0;%f2=figure;
sortedDrugNames = {};
for iMech = 1:nMech
    curInx = find(mechanismPrimaryInx==iMech);
    [~,drugInx] = sort(nPumpByDrug(curInx),'descend');
    junk = drugNames(curInx);
    sortedDrugNames = cat(1,sortedDrugNames,junk(drugInx));
    subplot(5,1,1:4); hold on;
    plot(lastN+1:lastN+numel(curInx),nPumps+1,'s','MarkerFaceColor',mechColor(iMech,:),'MarkerEdgeColor','none','MarkerSize',10)
    m=nPumps;
    for iPump = 1:nPumps
        curPump = sortedPumpImpact(iPump,curInx(drugInx));
        if (~isempty(find(curPump>0)))
            plot(find(curPump>0)+lastN,m,'o','MarkerFaceColor','w','MarkerEdgeColor','k','MarkerSize',5)
        end
        if (~isempty(find(curPump<0)))
            plot(find(curPump<0)+lastN,m,'o','MarkerFaceColor','k','MarkerEdgeColor','none','MarkerSize',5)
        end
        m=m-1;
    end
    subplot(5,1,5); hold on;
    curN = nPumpByDrug(curInx);
    bar(lastN+1:lastN+numel(curInx),curN(drugInx),'facecolor',mechColor(iMech,:),'BarWidth',1)
    lastN=lastN+(numel(curInx));
end

subplot(5,1,1:4);grid on;box on; ylim([0 nPumps+1]);yticks(1:nPumps+1);yticklabels(sortedPumpNames(nPumps:-1:1));%yticklabels([sortedPumpNames(nPumps:-1:1) 'mech' ]);
subplot(5,1,5);grid on;box on;P
set(gcf,'position',[1 484 1512 382]);
figure;bar(sum(sortedPumpImpact,2)*-1);grid on; box on;
% saveas(f1,'pumpsByMechSortedByfam.svg')
%% numbers
pumpsWithImpact = sum(sum(sortedPumpImpact,2)~=0)
sort(sum(sortedPumpImpact,1))
