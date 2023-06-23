function [tf xp yp] = getNodeInsidePolygon(fh,xq,yq,colorLabel)

if(nargin==3) % generate a random color
    myColors = jet;
    colorLabel = myColors(randi([1 256]),:);
end

[xp,yp] = ginput; % get points of polygon


plot([xp;xp(1)],[yp;yp(1)],'-','Color',colorLabel,'LineWidth',1.5);
pause(2);
tf = inpolygon(xq,yq,xp,yp);
plot(xq(tf),yq(tf),'o','MarkerFaceColor','none','MarkerEdgeColor',colorLabel,'MarkerSize',15,'LineWidth',1.5);

end