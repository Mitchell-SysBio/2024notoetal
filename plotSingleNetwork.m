function [fh,x,y,G,h,sim] = plotSingleNetwork(LFC,strainNames,drugNames,drugColors,strTitle,drugMarkerArray,corrCutoff,simMatrix,TFallEdges,TFforCytoscape)

if(nargin==6), corrCutoff = 0.45; end
if(nargin<9), TFallEdges=0; TFforCytoscape=0; end

% before anything, impute missing values by k nearest rows (strains)
data = LFC';
imputedData = knnimpute(data,5);
LFC2 = imputedData';

r = corr(LFC2,'rows','pairwise');
network = r;
if TFallEdges
else
    network(network<corrCutoff | network==1) = 0; % remove edges below correlation threshold and the diagonal (self-edges)
end
G = graph(network,'upper','omitselfloops');
fh = figure('color','white'); hold on;
h = plot(G,'NodeLabel',drugNames);
layout(h,'force','UseGravity',true,'weighteffect','inverse');
for i=1:length(network)
    highlight(h, i, 'NodeColor', drugColors(i,:),'MarkerSize',17);
end

h.Marker = drugMarkerArray;

sim = zeros(height(G.Edges),1);
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
        sim(i,1) = curSim;
        h.EdgeAlpha=1;
    end

end

h.LineWidth = ((G.Edges.Weight-min(G.Edges.Weight))./(max(G.Edges.Weight)-min(G.Edges.Weight))+0.1)*5;


x = h.XData; y = h.YData;
% add node numbers

if(isnumeric(drugNames))
    text(x,y,num2str(drugNames),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',8,'FontWeight','bold','Color','white');
end
% G.Nodes.Names = drugNames;

if TFforCytoscape;
    G = graph(network,drugNames,'upper','omitselfloops');
end;

title(strTitle);
set(gca,'xtick',[],'ytick',[]); box on; axis square; axis tight;
set(gcf,'Position',[230   208   806   658])
end