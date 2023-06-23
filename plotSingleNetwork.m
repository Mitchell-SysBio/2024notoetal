function [fh,x,y,G] = plotSingleNetwork(LFC,strainNames,drugNames,drugColors,strTitle,drugMarkerArray,corrCutoff,simMatrix)

if(nargin==6), corrCutoff = 0.45; end

% before anything, impute missing values by k nearest rows (strains)
data = LFC';
imputedData = knnimpute(data,5);
LFC2 = imputedData';

r = corr(LFC2,'rows','pairwise');
network = r;
network(network<corrCutoff | network==1) = 0; % remove edges below correlation threshold and the diagonal (self-edges)

G = graph(network,'upper');
fh = figure('color','white'); hold on;
h = plot(G,'NodeLabel','');
layout(h,'force','UseGravity',true,'weighteffect','inverse');
for i=1:length(network)
    highlight(h, i, 'NodeColor', drugColors(i,:),'MarkerSize',17);
end

h.Marker = drugMarkerArray;


if(nargin>6)
    junk = G.Edges;
    for i=1:height(junk)
        curJunk = junk(i,:);
        curSim = simMatrix(curJunk.EndNodes(1),curJunk.EndNodes(2));
        if curSim<0.2
            highlight(h, curJunk.EndNodes(1),curJunk.EndNodes(2), 'EdgeColor',[39,170,225]./255);
        else
            highlight(h, curJunk.EndNodes(1),curJunk.EndNodes(2), 'EdgeColor','k');       
        end
h.EdgeAlpha=1;
    end

end

h.LineWidth = ((G.Edges.Weight-min(G.Edges.Weight))./(max(G.Edges.Weight)-min(G.Edges.Weight))+0.1)*5;


x = h.XData; y = h.YData;

title(strTitle);
set(gca,'xtick',[],'ytick',[]); box on; axis square; axis tight;
set(gcf,'Position',[230   208   806   658])
end