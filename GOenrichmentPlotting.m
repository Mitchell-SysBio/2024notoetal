%% plots for figure 2B and 3C
% this script takes the gene ontology enrichment results for the drug
% groups in the antibiotic network or the non-antibiotic network. It makes
% an input for the online tool REVIGO and plots the output to show the top
% enriched pathways in each group
% 02-08-2023
%% load non-antibiotic clusters and pw annotation
load htClusters06082023.mat
goBpAnn = extractPwAnnFromGMT('./metadata/GO bp Entrez GMT.csv',5,1500); 
%% extract details of pw analysis 
% extract data of GO biological process enriched pathways
networkOrder = [1:5]; % clusters 6 & 7 are just subgroups of 3 & 5
enrichedPW = {}; enrichedCluster = []; enrichedPval = []; enrichedSet = {};enrichedUniv={};clusterSize=[];
for j = 1:numel(networkOrder);
    i = networkOrder(j);
    curPws = htClusters.enrichment.PwName{1,i}; 
    curPval = htClusters.enrichment.pvalFDR{1,i};
    curSet = htClusters.enrichment.genesSet{1,i};
    curUniv = htClusters.enrichment.genesUniverse{1,i};
    if length(curPws)>5
        enrichedPW = cat(1,enrichedPW,curPws(1:5));
        enrichedCluster = [enrichedCluster linspace(i,i,5)];
        enrichedPval = [enrichedPval; curPval(1:5)];
        enrichedSet = cat(1,enrichedSet,curSet(1:5));
        enrichedUniv = cat(1,enrichedUniv,curUniv(1:5));
        clusterSize = [clusterSize 5];
    else
        enrichedPW = cat(1,enrichedPW,curPws);
        enrichedCluster = [enrichedCluster linspace(i,i,length(curPws))];
        enrichedPval = [enrichedPval; curPval];
        enrichedSet = cat(1,enrichedSet,curSet);
        enrichedUniv = cat(1,enrichedUniv,curUniv);
        clusterSize = [clusterSize length(curPws)];
    end
end

%% make input for online Revigo analysis
% revigo compacts gene ontology results into general categories
clusterLabels = {'cellular biosynthetic processes' 'lipopolysaccharide biosynthesis' 'DNA damage' 'siderophore synthesis' 'oxidative stress'};
for iCluster = 1:5 % clusters with enriched processes
    inx = networkOrder(iCluster);
    curEnriched = htClusters.enrichment.PwName{1,inx};
    curPval = htClusters.enrichment.pvalFDR{1,inx};
    test={};
    for i = 1:length(curEnriched)
        curJunk = strsplit(curEnriched{i});
        test{i,1} = curJunk{1};
        test{i,2} = curPval(i);
    end
    writecell(test,'forRevigoHt0608.xlsx','Sheet',clusterLabels{iCluster});
end
%% do the same for antibiotics
load abClusters02092023.mat
% make revigo input 
for inx = 1:numel(abClusters.mechInx)
    curEnriched = abClusters.enrichment.PwName{1,inx};
    curPval = abClusters.enrichment.pvalFDR{1,inx};
    test={};
    for i = 1:length(curEnriched)
        curJunk = strsplit(curEnriched{i});
        test{i,1} = curJunk{1};
        test{i,2} = curPval(i);
    end
    writecell(test,'forRevigoAb0608.xlsx','Sheet',uniqMechanism{inx});
end
% extract data of GO biological process enriched pathways
enrichedPW = {}; enrichedCluster = []; enrichedPval = []; enrichedSet = {};enrichedUniv={};clusterSize=[];
for i = 1:numel(abClusters.mechInx);
%     i = networkOrder(j);
    curPws = abClusters.enrichment.PwName{1,i};
    curPval = abClusters.enrichment.pvalFDR{1,i};
    curSet = abClusters.enrichment.genesSet{1,i};
    curUniv = abClusters.enrichment.genesUniverse{1,i};
    if length(curPws)>5
        enrichedPW = cat(1,enrichedPW,curPws(1:5));
        enrichedCluster = [enrichedCluster linspace(i,i,5)];
        enrichedPval = [enrichedPval; curPval(1:5)];
        enrichedSet = cat(1,enrichedSet,curSet(1:5));
        enrichedUniv = cat(1,enrichedUniv,curUniv(1:5));
        clusterSize = [clusterSize 5];
    else
        enrichedPW = cat(1,enrichedPW,curPws);
        enrichedCluster = [enrichedCluster linspace(i,i,length(curPws))];
        enrichedPval = [enrichedPval; curPval];
        enrichedSet = cat(1,enrichedSet,curSet);
        enrichedUniv = cat(1,enrichedUniv,curUniv);
        clusterSize = [clusterSize length(curPws)];
    end
end

%% REVIGO PLOTS - plot by AB or HT
% plots for figure 2B and 3C
revigoDir = dir('RevigoOut*');
maxPwToPlot = 6;
for iFile = 1:length(revigoDir)
    curFile = revigoDir(iFile).name;
    clusterNames = sheetnames(curFile);
    legend = {};
     k = 1; m = numel(clusterNames); f=figure('color','white');hold on;
    for iCluster = 1:numel(clusterNames)
    pwNames = readcell(curFile,'Sheet',iCluster,'Range','B2:B35');
    pwPval = 10.^(readmatrix(curFile,'Sheet',iCluster,'Range','C2:C35'));
    pwSize = 10.^(readmatrix(curFile,'Sheet',iCluster,'Range','D2:D35'));
    inxDel = isnan(pwPval); pwNames(inxDel,:)=[];
    [~,inxToPlot] = sort(pwPval,'ascend');
    if length(pwNames)>=maxPwToPlot
        for iPW =1:maxPwToPlot
            curQval = pwPval(inxToPlot(iPW));
            curSetSize = pwSize(inxToPlot(iPW));
            curColor = min([-log10(curQval)/3 1]);
            plot(k,m,'o','markersize',log10(curSetSize)*10,'markerfacecolor',[curColor curColor curColor],'markeredgecolor','none');
            k = k+1;
        end
        legend = cat(1,legend,pwNames(inxToPlot(1:maxPwToPlot)));
    else
        for iPW =1:length(pwNames)
            curQval = pwPval(inxToPlot(iPW));
            curSetSize = pwSize(inxToPlot(iPW));
            curColor = min([-log10(curQval)/3 1]);
            plot(k,m,'o','markersize',log10(curSetSize)*10,'markerfacecolor',[curColor curColor curColor],'markeredgecolor','none');
            k = k+1;
        end
        legend = cat(1,legend,pwNames(inxToPlot(1:length(pwNames))));
    end
    m = m-1;
    end
    set(gca,'xtick',1:length(legend),'xticklabel',legend,'ytick',1:length(clusterNames),'yticklabel',clusterNames(numel(clusterNames):-1:1));
    set(gca,'xlim',[0.5 length(legend)+0.5],'ylim',[0 length(clusterNames)+1]);
    set(gcf,'position',[1 341 1512 457]);
    ytickangle(0); xtickangle(45);
    box on; grid on;
% saveas(f,['enrichmentPlots/' curFile(1:8) '.svg'])
% plot a "legend" image
fake_qval = [0.001,0.01,0.1 0.25];
fake_setSize = [10,50,250,500];
f = figure; hold on;
for i=1:4
curQval = fake_qval(i);
                curSetSize = fake_setSize(i);
                curColor = min([-log10(curQval)/3 1]);
                plot(i,1,'o','markersize',log10(curSetSize)*10,'markerfacecolor',[curColor curColor curColor],'markeredgecolor','none');
end
set(gca,'xtick',1:length(legend),'xticklabel',legend,'ytick',1:length(clusterNames),'yticklabel',clusterNames);
set(gca,'xlim',[0.5 length(legend)+0.5],'ylim',[0.5 length(clusterNames)+0.5]);
set(gcf,'position',[1 341 1512 457]);
ytickangle(0); xtickangle(45);
box on; grid on;
% saveas(f,['enrichmentPlots/' curFile(1:8) '_legend.svg'])
end
