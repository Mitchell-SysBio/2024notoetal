%% Load annotation data
load strainResultsAggregated.mat; %
raw = readcell('metadata/drugMetadata.xlsx'); %
%% Load screen data
% The raw drug metadata file is not in the same order as the screens in
% strainResults and includes ignored screens so we need to sort it first
% and then filter
screenID = strainResults.screenID'; % gets screen IDs
nScreens = length(screenID); % useful variable
screenInx2 = zeros(length(screenID),1);
for i=1:nScreens
    curID = screenID{i};
    %     curID = curID(1:regexp(curID,'_')-1);
    if ismember(i,[166:167])
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
drugMetadata = raw([screenInx2],:); % get medatadata on drugs sorted & filtered
drugMetadataHeader = raw(1,:);
%% Get drug category
drugNames = drugMetadata(:,strcmp(drugMetadataHeader,'Name')); % array to label with drug names
tfHt = logical(cell2mat(drugMetadata(:,strcmp(drugMetadataHeader,'antihost')))); % only host-targeted drugs
tfAbAll = logical(cell2mat(drugMetadata(:,strcmp(drugMetadataHeader,'antibiotic')))); % include unclassified antibiotics
drugCat = zeros(size(tfAbAll));
drugCat(tfAbAll) = 1;
drugCat(tfHt) = 2; % drug categories where 1=ab and 2=ht
%% load similarity data and sort
% similarity scores are assigned to drugs based on CID, this script matches
% the CID in the similarity score file with the CIDs in the metadata file
drugCID = [drugMetadata{:,contains(drugMetadataHeader,'CID')}]';
rawSim = readmatrix('drugNamesForQuery_simmat_upd100223_cid.xlsx');
inxCID = zeros(numel(drugCID),2);
for i = 1:numel(drugCID)
    curCID = drugCID(i);
    inxCID(i,1) = find(rawSim(:,1)==curCID);
    inxCID(i,2) = find(rawSim(1,:)==curCID);
end
similarityScores = rawSim(inxCID(:,1),inxCID(:,2));
%% separate between ab and ht
% makes an array with similarity scores of pairs between antibiotic and non-antibiotic
% drugs
abvsht=[];
for iRow = 1:nScreens
    curRowCat = drugCat(iRow);
    for iCol = iRow+1:nScreens
        curColCat = drugCat(iCol);
        if (curRowCat==curColCat)
        else
            abvsht = [abvsht, similarityScores(iRow,iCol)];
        end
    end
end
%% get mechanism for screens
mechTerm = 'mechanism'; % 'mechanism' for targeted mechanism
drugMechanismRaw = drugMetadata(:,strcmp(drugMetadataHeader,mechTerm));
mechanismPrimary = cell(nScreens,1); mechanismPrimaryInx = nan(nScreens,1);
uniqMechanism = unique(drugMechanismRaw(tfAbAll)); % unique mechanism terms
uniqMechanism{end+1} = 'anti-host';
mechanismPrimaryCounts = zeros(size(uniqMechanism)); % number of drugs per mechanism
for i = 1:nScreens
    curMech = drugMechanismRaw{i};
    if(tfHt(i)) % antihost drug
        mechanismPrimary{i} = 'anti-host';
    else % antibiotic drug
        mechanismPrimary{i} = drugMechanismRaw{i};
    end
    curPrimaryInx = find(strcmp(uniqMechanism,mechanismPrimary{i}));
    mechanismPrimaryInx(i) = curPrimaryInx; % for each drug what mechanism they're part of
    mechanismPrimaryCounts(curPrimaryInx) = mechanismPrimaryCounts(curPrimaryInx)+1;
end
mechColor = [186,85,211; 125,25,207; 220,20,60; 255,140,0; 255,215,0; 0,100,100; 0,128,0; 0,0,180; 240,128,128; 218,112,214; 100,100,100];
mechColor = mechColor./255;
mechColor(end,:) =  [0.7 0.7 0.7];
mechColor(end+1,:) =  [0.3 0.3 0.3]; % colormap for each mechanism
drugColorLabel = zeros(length(mechanismPrimaryInx),3); % colormap for each drug
for i=1:size(drugColorLabel,1)
    drugColorLabel(i,:) = mechColor(mechanismPrimaryInx(i),:);
end
clear drugMechanismRaw curPrimaryInx
%% plot histogram for fig 1E
toPlot=[]; % scores of antibiotic pairs with same mechanisms
btwGroup = []; % scores of antibiotic pairs with different mechanisms
for iRow = 1:nScreens
    curRowCat = drugCat(iRow);
    if curRowCat==1
        for iCol = iRow+1 :nScreens
            curColCat = drugCat(iCol);
            if (curRowCat==curColCat)
                rowMech = mechanismPrimaryInx(iRow);
                colMech = mechanismPrimaryInx(iCol);
                if (rowMech==colMech)
                    toPlot = [toPlot;similarityScores(iRow,iCol) ];
                else
                    btwGroup = [btwGroup;similarityScores(iRow,iCol) ];
                end
            else
                continue
            end
        end
    else
        continue
    end
end
f=figure('Color','white'); edges = [0:0.0333:1];
subplot(3,1,1);histogram(toPlot(:,1),edges,'Normalization','probability');grid on; box on;ylim([0 0.4])
subplot(3,1,2);histogram(btwGroup(:,1),edges,'Normalization','probability');grid on; box on;ylim([0 0.4])
subplot(3,1,3);histogram(abvsht,edges,'Normalization','probability');grid on; box on;ylim([0 0.4])
set(gcf,'position',[578 255 270 515])
%% cluster similarity scores and plot heatmap
d2 = pdist(similarityScores,'euclidean');
l2 = linkage(d2,'average');
inxToPlot = optimalleaforder(l2,d2);
markerColorToPlot = mechanismPrimaryInx(inxToPlot);
% make figure with clustered similarity scores and color labels for drugs
figure('Color','white');
junk=similarityScores(inxToPlot,inxToPlot);
junkNames = drugNames(inxToPlot);
imagesc(tril(junk),[0 1]);colormap(flipud(bone));
hold on;
for j = 1:numel(markerColorToPlot)
    plot(j,length(similarityScores(inxToPlot,inxToPlot))+0.5,'s','MarkerFaceColor',mechColor(markerColorToPlot(j),:),'MarkerEdgeColor',mechColor(markerColorToPlot(j),:),'MarkerSize',4)
    plot(0.5,j,'s','MarkerFaceColor',mechColor(markerColorToPlot(j),:),'MarkerEdgeColor',mechColor(markerColorToPlot(j),:),'MarkerSize',4)
end
axis square;
%% sort by mechanism AND within mechanism by score
[markerColorToPlot,inxToPlot] = sort(mechanismPrimaryInx,'ascend');
junk=similarityScores(inxToPlot,inxToPlot);
myStart = 1; medSim =[];
inxToPlot2 = zeros(size(inxToPlot));
for i=1:length(mechanismPrimaryCounts)
    curCounts = mechanismPrimaryCounts(i);
    myEnd = myStart-1+curCounts;
    curScores = junk(myStart:myEnd,myStart:myEnd);
    curInx = myStart:myEnd;
    if i==6
        lo2 = 1;
    else
        d2 = pdist(curScores,'euclidean');
        l2 = linkage(d2,'ward');
        lo2 = optimalleaforder(l2,d2);
    end
    inxToPlot2(myStart:myEnd) = curInx(lo2);
    myStart=myStart+curCounts;
    medSim = [medSim mean(curScores(:))];
end
figure('Color','white');
junk=junk(inxToPlot2,inxToPlot2);
imagesc(tril(junk),[0 1]);colormap(flipud(bone));
hold on;
for j = 1:numel(markerColorToPlot)
    plot(j,0.5,'s','MarkerFaceColor',mechColor(markerColorToPlot(j),:),'MarkerEdgeColor',mechColor(markerColorToPlot(j),:),'MarkerSize',4)
    plot(0.5,j,'s','MarkerFaceColor',mechColor(markerColorToPlot(j),:),'MarkerEdgeColor',mechColor(markerColorToPlot(j),:),'MarkerSize',4)
end
axis square;
% close all
%%
tfUnknown = logical(strcmp(drugMetadata(:,9),'unknown')); % drugs with unknown class/mechanism
tfAb = [tfAbAll & ~tfUnknown];
inxHt = find(tfHt);
toPlot=[]; % scores of antibiotic pairs with same mechanisms
btwGroup = []; % scores of antibiotic pairs with different mechanisms
for iRow = inxHt
    curRowCat = drugCat(iRow);
    curScores = similarityScores(iRow,tfAb);
    toPlot = [toPlot;max(curScores)];
end

abAb=zeros(nScreens,1); % scores of antibiotic pairs with same mechanisms
btwGroup = zeros(nScreens,1); % scores of antibiotic pairs with different mechanisms
for iRow = 1:nScreens
    curRowCat = drugCat(iRow);
    if curRowCat==1
        rowMech = mechanismPrimaryInx(iRow);
        curScores = similarityScores(iRow,mechanismPrimaryInx==rowMech);
        junk = sort(curScores,'descend');
        if numel(junk)>1
            abAb(iRow) = junk(2);
        end
        curScores = similarityScores(iRow,mechanismPrimaryInx~=rowMech);
        btwGroup(iRow) = max(curScores);
    end
end
abAb(abAb==0,:) = [];
btwGroup(btwGroup==0,:) = [];
f=figure('Color','white'); edges = [0:0.033:1];
subplot(3,1,1);histogram(abAb,edges,'Normalization','probability');grid on; box on;ylim([0 0.4])
subplot(3,1,2);histogram(btwGroup,edges,'Normalization','probability');grid on; box on;ylim([0 0.4])
subplot(3,1,3);histogram(toPlot,edges,'Normalization','probability');grid on; box on;ylim([0 0.4])
