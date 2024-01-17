%% user definitions
cutoffpadj = 0.25; % threshold to set  hits based on p-val
%% Load annotation data
load strainResultsAggregated.mat; % sequencing results for all screens
raw = readcell('metadata/drugMetadata.xlsx'); % drug info
strainMetadata = readcell('metadata/strainMetadata.csv'); % strain info
% e coli pathway annotations - convert entrez ID to names used in screen
goBpAnn = extractPwAnnFromGMT('./metadata/GO bp Entrez GMT.csv',5,1500); 
goCcAnn = extractPwAnnFromGMT('./metadata/GO cc Entrez GMT.csv',5,1500); 
goMfAnn = extractPwAnnFromGMT('./metadata/GO mf Entrez GMT.csv',5,1500); 
keggAnn = extractPwAnnFromGMT('./metadata/KEGG Entrez GMT.csv',5,1000); 
ecocycAnn = extractPwAnnFromGMT('./metadata/Ecocyc Pathways Entrez GMT.csv',5,1500);
cogAnn = extractPwAnnFromGMT('./metadata/COG function Entrez GMT.csv',5,1500);
%% Load screen data
% The raw drug metadata file is not in the same order as the screens in
% strainResults and includes ignored screens so we need to sort it first
% and then filter
screenID = strainResults.screenID'; % gets screen IDs
nScreens = length(screenID); % useful variable
screenInx2 = zeros(length(screenID),1); % index to match metadata to screen results
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
drugMetadataHeader = raw(1,:); % names of varibles in metadata

%% Get LFC and LFC masked by padj
LFC = strainResults.LFC; % LFC of filtered screens
padj = strainResults.padj; % padj of filtered screens
TFpadj = padj < cutoffpadj; % TF of screens with padj over threshold

%% standardize screens
stdLFC = nanstd(LFC,0,1);
normLFC = LFC./stdLFC;

% Parse annotation data
mechTerm = 'mechanism'; % 'mechanism' for targeted mechanism, 'class' for chemical class
inxMech = find(strcmp(drugMetadataHeader,mechTerm));
inxName = find(strcmp(drugMetadataHeader,'Name'));

% get TF where 1 means the screen was with an antibiotic, non-antibiotic or
% antibiotic with unknown mechanism of action 
tfAbAll = logical(cell2mat(drugMetadata(:,strcmp(drugMetadataHeader,'antibiotic')))); % include unclassified antibiotics
tfHt = logical(cell2mat(drugMetadata(:,strcmp(drugMetadataHeader,'antihost')))); % only host-targeted drugs
tfUnknown = logical(strcmp(drugMetadata(:,inxMech),'unknown')); % drugs with unknown class/mechanism
tfAb = [tfAbAll & ~tfUnknown];

drugNames = drugMetadata(:,inxName); % array to label with drug names
drugMechanismRaw = drugMetadata(:,inxMech);
uniqMechanism = unique(drugMechanismRaw(tfAbAll)); % unique mechanism terms
uniqMechanism{end+1} = 'anti-host'; 
%% Parse mechanisms and prepare mechanism arrays
mechanismPrimary = cell(nScreens,1); mechanismPrimaryInx = nan(nScreens,1); 
mechanismPrimaryCounts = zeros(size(uniqMechanism)); % number of drugs per mechanism

for i = 1:nScreens
    curMech = drugMechanismRaw{i};
    if(tfHt(i)) % antihost drug
        mechanismPrimary{i} = 'anti-host'; 
    else % antibiotic drug
        mechanismPrimary{i} = drugMechanismRaw{i};
    end
    curPrimaryInx = find(strcmp(uniqMechanism,mechanismPrimary{i}));
    mechanismPrimaryInx(i) = curPrimaryInx;
    mechanismPrimaryCounts(curPrimaryInx) = mechanismPrimaryCounts(curPrimaryInx)+1;
end

% house cleaning
clear raw i curID inxMech inxName uniqAbMechanism curPrimaryInx screenInx2;
%% Parse non-antibiotic mechanisms and prepare mechanism arrays
% drugMechanismRaw = drugMetadata(:,strcmp(drugMetadataHeader,mechTerm));
uniqMechanismHt = unique(drugMechanismRaw(tfHt)); % unique mechanism terms
uniqMechanismHt{end+1} = 'ab'; 
uniqMechanismHt{end+1} = 'ab-unknown'; 

mechanismHt = cell(nScreens,1); mechanismHtInx = nan(nScreens,1); 
mechanismHtCounts = zeros(size(uniqMechanismHt)); % number of drugs per mechanism
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

%% Build colormap
% colormap for antibiotic mechanisms and black for non-antibiotic
mechColor = [186,85,211; 125,25,207; 220,20,60; 255,140,0; 255,215,0; 0,100,100; 0,128,0; 0,0,180; 240,128,128; 218,112,214; 100,100,100];
mechColor = mechColor./255;
mechColor(end,:) =  [0.7 0.7 0.7]; % unknown antibiotic
mechColor(end+1,:) =  [0.3 0.3 0.3]; % non-antibiotic
drugColorLabel = zeros(length(mechanismPrimaryInx),3);
for i=1:size(drugColorLabel,1)
    drugColorLabel(i,:) = mechColor(mechanismPrimaryInx(i),:);
end

fh=figure('color','white'); 
pie(mechanismPrimaryCounts,[zeros(1,11) 1],num2str(mechanismPrimaryCounts)); 
legend(uniqMechanism,'location','westoutside');
colormap(mechColor); 

% markers for non-antibiotic
markerArray = {'o' 'h' 'diamond' 'v' 's' '^' 'o' 'o'};
drugMarkerArray = {};
for i=1:nScreens
    drugMarkerArray(i,:) = markerArray(mechanismHtInx(i));
end

%% Infer useful indexes for drug groups and gene groups
% for filtering of strains
% find slow growers
growthRate= cell2mat(strainMetadata(ismember(strainMetadata(:,1),strainNames),17));
growthRateNames = strainMetadata(ismember(strainMetadata(:,1),strainNames),1);
slowGrowerNames = growthRateNames(growthRate<0.7);
inxSlowStrains = find(ismember(strainNames,slowGrowerNames));
% find prophage strains
load ./metadata/prophageStrains.mat
prophages(strcmp(prophages(:,1),'yagN'),:)=[];
inxProphageStrains = find(ismember(strainNames,prophages(:,1)));
% find eflux systems and membrane related strains
load ./metadata/pumpStrains.mat
inxPumpStrains = find(ismember(strainNames,pumps(:,1)));
inxMembraneStrains = unique([find(startsWith(strainNames,'omp'));find(startsWith(strainNames,'tol'));...
    find(startsWith(strainNames,'ton'));find(startsWith(strainNames,'pal'));find(startsWith(strainNames,'smpA'))...
    ;find(startsWith(strainNames,'yfgL'));find(startsWith(strainNames,'asmA'));find(startsWith(strainNames,'surA'));...
    ;find(startsWith(strainNames,'exbB'));find(startsWith(strainNames,'exbD'));find(startsWith(strainNames,'marR'));...
    find(startsWith(strainNames,'rob'));find(startsWith(strainNames,'hlpA')); find(startsWith(strainNames,'yrb'))]);%; find(startsWith(strainNames,'yag'))]);

inxAllDrugs = [1:nScreens]';
inxAbDrugs = find(tfAb);
inxHtDrugs = find(tfHt);
inxAbAllDrugs = find(tfAbAll);
inxAllCharDrugs = [inxAbDrugs; inxHtDrugs];

inxAllStrains = 1:size(normLFC,1);
inxRemoveStrains = [inxSlowStrains;inxProphageStrains;inxPumpStrains;inxMembraneStrains];
inxAllUsefullStrains = setdiff(inxAllStrains,inxRemoveStrains);

TFLFC = abs(normLFC) > log2(2.5);
tfHits = TFLFC & TFpadj;

%% load similarity data and sort
% added 02-14-2023
% similarity scores are assigned to drugs based on CID, this script matches
% the CID in the similarity score file with the CIDs in the drug metadata file
drugCID = [drugMetadata{:,contains(drugMetadataHeader,'CID')}]';
rawSim = readmatrix('drugNamesForQuery_simmat_upd100223_cid.xlsx');
inxCID = zeros(numel(drugCID),2);
for i = 1:numel(drugCID)
    curCID = drugCID(i);
    inxCID(i,1) = find(rawSim(:,1)==curCID);
    inxCID(i,2) = find(rawSim(1,:)==curCID);
end
similarityScores = rawSim(inxCID(:,1),inxCID(:,2));
%% Highly similar drugs
data = normLFC';
imputedData = knnimpute(data,5);
LFC2 = imputedData';
r = corr(LFC2,'rows','pairwise');
abColorLabel = drugColorLabel(inxAbDrugs,:);
simabht = similarityScores(inxAbDrugs,inxHtDrugs);
rabht = r(inxAbDrugs,inxHtDrugs);
f=figure('color','white'); hold on; plot(simabht(simabht<0.2),rabht(simabht<0.2),'o',...
    'MarkerEdgeColor',[0.7 0.7 0.7],'MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',3);
keepSim = simabht; keepSim(keepSim<0.2)=0;
k=1;
for i = 1:size(keepSim,1)
    for j = 1:size(keepSim,2)
        if keepSim(i,j)>0
            plot(simabht(i,j),rabht(i,j),'o','MarkerEdgeColor',abColorLabel(i,:),'MarkerFaceColor',abColorLabel(i,:),'MarkerSize',5)
k=k+1
        end 
    end
end
xlim([0 1]);ylim([-0.2 1]); axis square;
grid on; box on; ylabel('drug screen correlation'); xlabel('similarity score')
title('antibiotic - non-antibiotic pair');%legend({'below 0.2 threshold','above 0.2 threshold'})
saveas(f,[date 'chemicalSimilarity_vs_screenCorrelationABHT.svg'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot networks and hierarchical clustering %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inxAllDrugs = [1:nScreens]';
inxAbDrugs = find(tfAb);
inxHtDrugs = find(tfHt);
inxAbAllDrugs = find(tfAbAll);
inxAllCharDrugs = [inxAbDrugs; inxHtDrugs];

inxAllStrains = 1:size(normLFC,1);
inxRemoveStrains = [inxSlowStrains;inxProphageStrains;inxPumpStrains;inxMembraneStrains];
inxAllUsefullStrains = setdiff(inxAllStrains,inxRemoveStrains);

TFLFC = abs(normLFC) > log2(2.5);
tfHits = TFLFC & TFpadj;

%% Ab based analysis
nStrainHits = sum(tfHits(:,inxAbDrugs),2);
inxAbInformativeStrains = find(nStrainHits>3 & nStrainHits <30);
inxStrains = setdiff(inxAbInformativeStrains,inxRemoveStrains);
[fh,~,~,G] = plotSingleNetwork(normLFC(inxStrains,inxAbDrugs),strainNames(inxStrains),inxAbDrugs,drugColorLabel(inxAbDrugs,:),'Ab by Ab strains',drugMarkerArray(inxAbDrugs,:),0.45,similarityScores(inxAbDrugs,inxAbDrugs));
fh = plotSingleHCluster(normLFC(inxStrains,inxAbDrugs),strainNames(inxStrains),drugNames(inxAbDrugs),drugColorLabel(inxAbDrugs,:),'Ab by Ab strains',true);
[fh,~,~,G] = plotSingleNetwork(normLFC(inxStrains,inxAllDrugs),strainNames(inxStrains),inxAllDrugs,drugColorLabel(inxAllDrugs,:),'All drugs by Ab strains',drugMarkerArray(inxAllDrugs,:),0.45,similarityScores(inxAllDrugs,inxAllDrugs));
%% Cytoscape networks
% for cytoscape antibiotics with filtered interactions
[fh,~,~,G,h,sim] = plotSingleNetwork(normLFC(inxStrains,inxAbDrugs),strainNames(inxStrains),drugNames(inxAbDrugs),drugColorLabel(inxAbDrugs,:),'Ab by Ab strains',drugMarkerArray(inxAbDrugs,:),0.45,similarityScores(inxAbDrugs,inxAbDrugs),0,1);
cyto = G.Edges; 
edgeWeight = G.Edges.Weight; %cyto.Weight = edgeWeight';junk = h.EdgeColor; colorEdge = sum(junk,2)>0; cyto.SimColor = colorEdge;
cyto.SimColor = sim;
missingDrugs=setdiff(drugNames(inxAbDrugs),unique(cyto.EndNodes)); missingDrugs = cat(2,missingDrugs,missingDrugs); junk=linspace(0.01,0.1,length(missingDrugs));
junk2 = linspace(0,0,length(missingDrugs));
missingDrugs=table(missingDrugs,junk',junk2','VariableNames',{'EndNodes','Weight','SimColor'});
atributes = table(drugNames(inxAbDrugs),mechanismPrimary(inxAbDrugs),'VariableNames',{'drug', 'mechanism'});writetable(atributes,[date 'nodeCytoscapeAb.xlsx']);
cytoFinal = cat(1,cyto,missingDrugs); writetable(cytoFinal,[date 'edgesCytoscapeAb.xlsx']); 
% for cytoscape antibiotics with all interactions
[fh,~,~,G,h,sim] = plotSingleNetwork(normLFC(inxStrains,inxAbDrugs),strainNames(inxStrains),drugNames(inxAbDrugs),drugColorLabel(inxAbDrugs,:),'Ab by Ab strains',drugMarkerArray(inxAbDrugs,:),0.45,similarityScores(inxAbDrugs,inxAbDrugs),1,1);
cyto = G.Edges; 
edgeWeight = G.Edges.Weight; %cyto.Weight = edgeWeight'; junk = h.EdgeColor; colorEdge = sum(junk,2)>0; cyto.SimColor = colorEdge;
cyto.SimColor = sim;
missingDrugs=setdiff(drugNames(inxAbDrugs),unique(cyto.EndNodes)); missingDrugs = cat(2,missingDrugs,missingDrugs); junk=linspace(0.01,0.1,length(missingDrugs));
junk2 = linspace(0,0,length(missingDrugs));
missingDrugs=table(missingDrugs,junk',junk2','VariableNames',{'EndNodes','Weight','SimColor'});
atributes = table(drugNames(inxAbDrugs),mechanismPrimary(inxAbDrugs),'VariableNames',{'drug', 'mechanism'});writetable(atributes,[date 'nodeCytoscapeAbAllEdges.xlsx']);
cytoFinal = cat(1,cyto,missingDrugs); writetable(cytoFinal,[date 'edgesCytoscapeAbAllEdges.xlsx']); 
% all drugs based on antibiotic analysis with filtered interactions
[fh,~,~,G,h,sim] = plotSingleNetwork(normLFC(inxStrains,inxAllDrugs),strainNames(inxStrains),drugNames(inxAllDrugs),drugColorLabel(inxAllDrugs,:),'All drugs by Ab strains',drugMarkerArray(inxAllDrugs,:),0.45,similarityScores(inxAllDrugs,inxAllDrugs),0,1);
cyto = G.Edges; 
edgeWeight = G.Edges.Weight; %cyto.Weight = edgeWeight'; junk = h.EdgeColor; colorEdge = sum(junk,2)>0; cyto.SimColor = colorEdge;
cyto.SimColor = sim;
missingDrugs=setdiff(drugNames(inxAllDrugs),unique(cyto.EndNodes)); missingDrugs = cat(2,missingDrugs,missingDrugs); junk=linspace(0.01,0.1,length(missingDrugs));
junk2 = linspace(0,0,length(missingDrugs));
missingDrugs=table(missingDrugs,junk',junk2','VariableNames',{'EndNodes','Weight','SimColor'});
atributes = table(drugNames(inxAllDrugs),mechanismPrimary(inxAllDrugs),'VariableNames',{'drug', 'mechanism'});writetable(atributes,[date 'nodeCytoscapeAbAlldrugs.xlsx']);
cytoFinal = cat(1,cyto,missingDrugs); writetable(cytoFinal,[date 'edgesCytoscapeAbAllDrugs.xlsx']); 
% all drugs based on antibiotic analysis with all interactions
[fh,~,~,G,h,sim] = plotSingleNetwork(normLFC(inxStrains,inxAllDrugs),strainNames(inxStrains),drugNames(inxAllDrugs),drugColorLabel(inxAllDrugs,:),'All drugs by Ab strains',drugMarkerArray(inxAllDrugs,:),0.45,similarityScores(inxAllDrugs,inxAllDrugs),1,1);
cyto = G.Edges; 
edgeWeight = G.Edges.Weight; %cyto.Weight = edgeWeight'; junk = h.EdgeColor; colorEdge = sum(junk,2)>0; cyto.SimColor = colorEdge;
cyto.SimColor = sim;
missingDrugs=setdiff(drugNames(inxAllDrugs),unique(cyto.EndNodes)); missingDrugs = cat(2,missingDrugs,missingDrugs); junk=linspace(0.01,0.1,length(missingDrugs));
junk2 = linspace(0,0,length(missingDrugs));
missingDrugs=table(missingDrugs,junk',junk2','VariableNames',{'EndNodes','Weight','SimColor'});
atributes = table(drugNames(inxAllDrugs),mechanismPrimary(inxAllDrugs),'VariableNames',{'drug', 'mechanism'});writetable(atributes,[date 'nodeCytoscapeAbAlldrugsAllEdges.xlsx']);
cytoFinal = cat(1,cyto,missingDrugs); writetable(cytoFinal,[date 'edgesCytoscapeAbAllDrugsAllEdges.xlsx']); 
%% Ht based analysis
nStrainHits = sum(tfHits(:,inxHtDrugs),2);
inxHtInformativeStrains = find(nStrainHits>2 & nStrainHits <30);
inxStrains = setdiff(inxHtInformativeStrains,inxRemoveStrains);
fh = plotSingleNetwork(normLFC(inxStrains,inxHtDrugs),strainNames(inxStrains),inxHtDrugs,drugColorLabel(inxHtDrugs,:),'Ht by Ht strains',drugMarkerArray(inxHtDrugs,:),0.46,similarityScores(inxHtDrugs,inxHtDrugs));
% saveas(fh,'Ht network by Ht strains.svg')
fh = plotSingleHCluster(normLFC(inxStrains,inxHtDrugs),strainNames(inxStrains),inxHtDrugs,drugColorLabel(inxHtDrugs,:),'Ht by Ht strains',true);
%% Cytoscape networks
%for cytoscape non-antibiotics with filtered interactions
[fh,~,~,G,h,sim] = plotSingleNetwork(normLFC(inxStrains,inxHtDrugs),strainNames(inxStrains),drugNames(inxHtDrugs),drugColorLabel(inxHtDrugs,:),'Ab by Ab strains',drugMarkerArray(inxHtDrugs,:),0.46,similarityScores(inxHtDrugs,inxHtDrugs),0,1);
cyto = G.Edges; 
edgeWeight = G.Edges.Weight; %cyto.Weight = edgeWeight';junk = h.EdgeColor; colorEdge = sum(junk,2)>0; cyto.SimColor = colorEdge;
cyto.SimColor = sim;
missingDrugs=setdiff(drugNames(inxHtDrugs),unique(cyto.EndNodes)); missingDrugs = cat(2,missingDrugs,missingDrugs); junk=linspace(0.01,0.1,length(missingDrugs));
junk2 = linspace(0,0,length(missingDrugs));
missingDrugs=table(missingDrugs,junk',junk2','VariableNames',{'EndNodes','Weight','SimColor'});
atributes = table(drugNames(inxHtDrugs),mechanismHt(inxHtDrugs),'VariableNames',{'drug', 'mechanism'});
writetable(atributes,[date 'nodeCytoscapeHt.xlsx']);
cytoFinal = cat(1,cyto,missingDrugs); writetable(cytoFinal,[date 'edgesCytoscapeHt.xlsx']); 
%for cytoscape non-antibiotics with all interactions
[fh,~,~,G,h,sim] = plotSingleNetwork(normLFC(inxStrains,inxHtDrugs),strainNames(inxStrains),drugNames(inxHtDrugs),drugColorLabel(inxHtDrugs,:),'Ab by Ab strains',drugMarkerArray(inxHtDrugs,:),0.46,similarityScores(inxHtDrugs,inxHtDrugs),1,1);
cyto = G.Edges; 
edgeWeight = G.Edges.Weight; %cyto.Weight = edgeWeight';junk = h.EdgeColor; colorEdge = sum(junk,2)>0; cyto.SimColor = colorEdge;
cyto.SimColor = sim;
missingDrugs=setdiff(drugNames(inxHtDrugs),unique(cyto.EndNodes)); missingDrugs = cat(2,missingDrugs,missingDrugs); junk=linspace(0.01,0.1,length(missingDrugs));
junk2 = linspace(0,0,length(missingDrugs));
missingDrugs=table(missingDrugs,junk',junk2','VariableNames',{'EndNodes','Weight','SimColor'});
atributes = table(drugNames(inxHtDrugs),mechanismHt(inxHtDrugs),'VariableNames',{'drug', 'mechanism'});
writetable(atributes,[date 'nodeCytoscapeHtAll.xlsx']);
cytoFinal = cat(1,cyto,missingDrugs); writetable(cytoFinal,[date 'edgesCytoscapeHtAll.xlsx']); 
%% Ab based analysis - manually curated clusters (polygon) from network
nStrainHits = sum(tfHits(:,inxAbDrugs),2);
inxAbInformativeStrains = find(nStrainHits>3 & nStrainHits <30);
inxStrains = setdiff(inxAbInformativeStrains,inxRemoveStrains);

[fh x y] = plotSingleNetwork(normLFC(inxStrains,inxAbDrugs),strainNames(inxStrains),inxAbDrugs,drugColorLabel(inxAbDrugs,:),'Ab by Ab strains',drugMarkerArray(inxAbDrugs));
set(gcf,'position',[260   328   907   619]);
abClusters.mechInx = [1 2 3 4 5 7 9];
abClusters.mechRelevant = uniqMechanism(abClusters.mechInx);
drugNamesAb = drugNames(inxAbDrugs);
for i=1:numel(abClusters.mechRelevant)
    disp(['Tag a polygon to define the ' abClusters.mechRelevant{i} ' cluster']);
    [abClusters.tfClusterDrugs{i} abClusters.xp{i} abClusters.yp{i}]  = getNodeInsidePolygon(fh,x,y,mechColor(abClusters.mechInx(i),:));
    abClusters.drugNames{i} = drugNamesAb(abClusters.tfClusterDrugs{i});
end
% saveas(fh,'Ab network by Ab strains - polygon.svg')
%% Plot subgraph for each cluster
normLFCab = normLFC(inxStrains,inxAbDrugs);
data = normLFCab';
imputedData = knnimpute(data,5);
normLFCab= imputedData';

for i=1:numel(abClusters.mechRelevant)
    curLFC = normLFCab(:,abClusters.tfClusterDrugs{i});
    curR = corr(curLFC);
    curTF = curR>0.45; % same cutoff used when building the network (plotSingleNetwork.m)
    curTF(1:size(curTF,2)+1:end)=0;
    G = graph(curTF,'upper');
    fh=figure('color','white');p=plot(G,'Layout','circle','NodeLabel','');
    curClusterName = abClusters.mechRelevant{i};
    p.NodeColor = mechColor(startsWith(uniqMechanism,curClusterName(1:3)),:);
    p.EdgeColor = 'k';
    title(curClusterName); axis square;
    % saveas(fh,[curClusterName '.svg'])
end
%% pathway enrichment analysis on manually curated clusters
pvalFDR_cutoff_strain = 0.1; % cutoff to define a significant enriched strains (from previous analysis)
pvalFDR_cutoff_pw = 0.25; % cutoff to define significantly enriched pathways
pwSizeRange = [10 1500];
nSet = length(inxStrains)/numel([abClusters.mechInx]);
ontologyArray = {goBpAnn;goCcAnn;goMfAnn;keggAnn;cogAnn;ecocycAnn};
ontologyDesc = {'goBpAnn'; 'goCcAnn'; 'goMfAnn'; 'keggAnn';'cogAnn';'ecocycAnn'};
for i=1:numel(abClusters.mechRelevant)
    [~,abClusters.pvalFDR{i}] = findHitStrains(normLFC(inxStrains,inxAbDrugs),inxAbDrugs,abClusters.tfClusterDrugs{i},pvalFDR_cutoff_strain);
    [~,inx] = sort(abClusters.pvalFDR{i});

    [~,inx] = sort(abClusters.pvalFDR{i},'ascend'); % add something that makes sure they're significant
    inxTF = inx(abClusters.pvalFDR{i}(inx)<pvalFDR_cutoff_strain);
    if length(inxTF)>nSet
            geneSet = strainNames(inxStrains(inxTF(1:nSet)));
    else
         geneSet = strainNames(inxStrains(inxTF));
    end
    geneUniverse = strainNames(inxStrains);
    % enrichment by goBp
    for j = 1:length(ontologyArray)
        if (j==6 | j==4)
            pwSizeRange = [5 50];
        end
        [PwName pvalFDR PwGenes PwGenesUniv PwGeneSet] = checkPwEnrichment(geneSet,geneUniverse,ontologyArray{j},pwSizeRange);
        inx = pvalFDR<pvalFDR_cutoff_pw;
        PwName = PwName(inx); pvalFDR = pvalFDR(inx); %overwrite variables keeping only significant ones
        [~,inx] = sort(pvalFDR,'ascend');
        abClusters.enrichment.pvalFDR{j,i} = pvalFDR(inx);
        abClusters.enrichment.PwName{j,i} = PwName(inx);
        abClusters.enrichment.OntologyName(j,1) = ontologyDesc(j);
        abClusters.enrichment.nSet(i) = numel(geneSet);
             abClusters.enrichment.nSet(i) = numel(geneSet);
       abClusters.enrichment.genesPresent{j,i} = PwGenes(inx)';
        abClusters.enrichment.genesUniverse{j,i} = PwGenesUniv(inx)';
        abClusters.enrichment.genesSet{j,i} = PwGeneSet(inx)';


    end

end
save abClusters06082023.mat abClusters
%% HT based analysis
TFLFC = abs(normLFC) > log2(2.5);
tfHits = TFLFC & TFpadj;
nStrainHits = sum(tfHits(:,inxHtDrugs),2);
inxHtInformativeStrains = find(nStrainHits>2 & nStrainHits <30);
inxStrains = setdiff(inxHtInformativeStrains,inxRemoveStrains);

[fh x y] = plotSingleNetwork(normLFC(inxStrains,inxHtDrugs),strainNames(inxStrains),inxHtDrugs,drugColorLabel(inxHtDrugs,:),'Ht by Ht strains',drugMarkerArray(inxHtDrugs),0.46,similarityScores(inxHtDrugs,inxHtDrugs));
set(gcf,'position',[260   328   907   619]);
nClusters = 9;
htClusters.mechInx = nan(nClusters,1);
htClusters.mechRelevant = nan(nClusters,1);
myJet = jet(nClusters);
drugNamesHt = drugNames(inxHtDrugs);
for i=1:nClusters
    disp(['Tag a polygon to define the ' num2str(i) ' cluster']);
    [htClusters.tfClusterDrugs{i} htClusters.xp{i} htClusters.yp{i}]  = getNodeInsidePolygon(fh,x,y,myJet(i,:));
    htClusters.drugNames{i} = drugNamesHt(htClusters.tfClusterDrugs{i});
end
% saveas(fh,'Ht network by Ht strains - polygon.svg')

%% Plot subgraph for each cluster
clusterNames = readcell('nodesTable.xlsx','Range','A2:A15');
normLFCht = normLFC(inxStrains,inxHtDrugs);
data = normLFCht';
imputedData = knnimpute(data,5);
normLFCht= imputedData';

for i=1:numel(htClusters.tfClusterDrugs)
    curLFC = normLFCht(:,htClusters.tfClusterDrugs{i});
    curR = corr(curLFC);
    curTF = curR>0.45; % same cutoff used when building the network (plotSingleNetwork.m)
    curTF(1:size(curTF,2)+1:end)=0;
    G = graph(curTF,'upper');
    fh=figure('color','white');p=plot(G,'Layout','circle','NodeLabel','');
    p.NodeColor = 'k';
    p.EdgeColor = 'k';
    curClusterName = clusterNames{i+7};
    title(curClusterName); axis square;
    % saveas(fh,[curClusterName '.svg'])
end
%%
pvalFDR_cutoff_strain = 0.1; % cutoff to define a significant enriched strains (from previous analysis)
pvalFDR_cutoff_pw = 0.25; % cutoff to define significantly enriched pathways
pwSizeRange = [10 1500];
nSet = length(inxStrains)/numel([htClusters.mechInx]);
ontologyArray = {goBpAnn;goCcAnn;goMfAnn;keggAnn;cogAnn;ecocycAnn};
ontologyDesc = {'goBpAnn'; 'goCcAnn'; 'goMfAnn'; 'keggAnn';'cogAnn';'ecocycAnn'};
for i=1:numel(htClusters.mechRelevant)
    [~,htClusters.pvalFDR{i}] = findHitStrains(normLFC(inxStrains,inxHtDrugs),inxHtDrugs,htClusters.tfClusterDrugs{i},pvalFDR_cutoff_strain);
    [~,inx] = sort(htClusters.pvalFDR{i});

    [~,inx] = sort(htClusters.pvalFDR{i},'ascend'); % add something that makes sure they're significant
    inxTF = inx(htClusters.pvalFDR{i}(inx)<pvalFDR_cutoff_strain);
    if length(inxTF)>nSet
            geneSet = strainNames(inxStrains(inxTF(1:nSet)));
    else
         geneSet = strainNames(inxStrains(inxTF));
    end
    geneUniverse = strainNames(inxStrains);
    % enrichment by goBp
    for j = 1:length(ontologyArray)
        if (j==6 | j==4)
            pwSizeRange = [5 50];
        end
        [PwName pvalFDR PwGenes PwGenesUniv PwGeneSet] = checkPwEnrichment(geneSet,geneUniverse,ontologyArray{j},pwSizeRange); %added 02-08-2023
        inx = pvalFDR<pvalFDR_cutoff_pw;
        PwName = PwName(inx); pvalFDR = pvalFDR(inx); %overwrite varihtles keeping only significant ones
        [~,inx] = sort(pvalFDR,'ascend');
        htClusters.enrichment.pvalFDR{j,i} = pvalFDR(inx);
        htClusters.enrichment.PwName{j,i} = PwName(inx);
        htClusters.enrichment.OntologyName(j,1) = ontologyDesc(j);
        htClusters.enrichment.nSet(i) = numel(geneSet);
        htClusters.enrichment.genesPresent{j,i} = PwGenes(inx)';
        htClusters.enrichment.genesUniverse{j,i} = PwGenesUniv(inx)';
        htClusters.enrichment.genesSet{j,i} = PwGeneSet(inx)';

    end

end
%%
save htClusters06082023.mat htClusters

save networkAnalysisVariables.mat drugNames inxAbDrugs inxHtDrugs normLFC tfHits inxRemoveStrains ...
    uniqMechanism mechColor mechanismPrimary drugColorLabel
save dataForRF.mat tfAb mechanismPrimary normLFC drugMechanismRaw inxAbAllDrugs...
tfAbAll mechColor tfHt inxHtDrugs strainNames drugColorLabel mechanismPrimaryInx

