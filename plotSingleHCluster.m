function [fh,lo2,lo1] = plotSingleHCluster(LFC,strainNames,drugNames,drugColors,strTitle,tfDendrogram)

if(nargin==5)
    tfDendrogram = false;
end

clim = [-3 3]; % color limits

% before anything, impute missing values by k nearest rows (strains)
if(sum(sum(isnan(LFC))))
    data = LFC';
    imputedData = knnimpute(data,5);
    LFC2 = imputedData';
else
    LFC2 = LFC;
end

d1 = pdist(LFC2,'euclidean');
d2 = pdist(LFC2','euclidean');
if (size(LFC,1)==1)
    l1 = 1;
    lo1 = l1;
else
    l1 = linkage(d1,'ward');
    lo1 = optimalleaforder(l1,d1);
end

l2 = linkage(d2,'ward');

lo2 = optimalleaforder(l2,d2);

fh = figure('color','white'); hold on;

if(tfDendrogram)
    subplot(6,1,1);
    dh = dendrogram(l2,0,'reorder',lo2);
    set(dh(:), 'Color', 'k' );
    set(gca,'XTick',[],'YTick',[],'Visible','off','box','off');
    subplot(6,1,[2:6]);
end

hm = imagesc(LFC2(lo1,lo2),clim);colormap(redbluecmap);
set(gca,'XTick',[1:numel(lo2)],'YTick',[1:numel(lo1)],'XTickLabel',drugNames(lo2),'YTickLabel',strainNames(lo1),'XTickLabelRotation',45)
hold on;
for j = 1:numel(lo2)
    plot(j,0.5,'s','MarkerFaceColor',drugColors(lo2(j),:),'MarkerEdgeColor',drugColors(lo2(j),:),'MarkerSize',10)
end
title(strTitle)
axis tight;

set(gcf,'position',[1          87        1512         779]);
end