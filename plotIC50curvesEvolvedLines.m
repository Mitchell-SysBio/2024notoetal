%% Plot IC50 curves evolved lines to Streptozotocin sertraline and triclabendazole
% 2023-03-02 MN
% fig 5b
load CrossResistanceABX.mat % output of Carmen's CrossResABXIC50Analysis384Well.m
%% concatenate all data
% Carmen saved data cells based on the plate things came from so I cat
% everything into a cell where each row is an individual drug-line pair
IC50data = {};
for i =1:length(IC50total)
    IC50data = cat(1,IC50data,table2cell(IC50total{i, 1}{1, 1}));
end

% junk = IC50data(contains(IC50data(:,5),'Sert') & contains(IC50data(:,1),'Sert'),:);
IC50data{102, 4}{1, 1}{1, 2} = [1;1.3411;0.9289;1.2682;0.4673;0.0297];
IC50data{102, 4}{2, 1}{1, 2} = [1;1.2456;0.8289;1.2216;0.5234;0.0355];
IC50data{102, 4}{3, 1}{1, 2} = [1;1.23140123804164;0.925492402926280;1.31716375914463;0.425773776027012;0.0280247608328644];
IC50data{115, 4}{3, 1}{1, 2}(end) = 0;IC50data{115, 4}{2, 1}{1, 2}(end) = 0;IC50data{115, 4}{1, 1}{1, 2}(end) = 0;
dose = [IC50data{102, 4}{1, 1}{1, 1},IC50data{102, 4}{2, 1}{1, 1},IC50data{102, 4}{3, 1}{1, 1}];
response = [IC50data{102, 4}{1, 1}{1, 2},IC50data{102, 4}{2, 1}{1, 2},IC50data{102, 4}{3, 1}{1, 2}];
[~,~,~,~,~, xpoints, ypoints] = calcDoseResponse(dose(:),response(:),false ,50);
IC50data{102, 7}{2, 1}{1, 1}  = xpoints;
IC50data{102, 7}{2, 1}{1, 2}  = ypoints;
IC50data{102, 7}{1, 1}{1, 1}  = xpoints;
IC50data{102, 7}{1, 1}{1, 2}  = ypoints;
IC50data{102, 7}{3, 1}{1, 1}  = xpoints;
IC50data{102, 7}{3, 1}{1, 2}  = ypoints;
dose = [IC50data{115, 4}{1, 1}{1, 1},IC50data{115, 4}{2, 1}{1, 1},IC50data{115, 4}{3, 1}{1, 1}];
response = [IC50data{115, 4}{1, 1}{1, 2},IC50data{115, 4}{2, 1}{1, 2},IC50data{115, 4}{3, 1}{1, 2}];
[~,~,~,~,~, xpoints, ypoints] = calcDoseResponse(dose(:),response(:),false ,50);
IC50data{115, 7}{2, 1}{1, 1}  = xpoints;
IC50data{115, 7}{2, 1}{1, 2}  = ypoints;
IC50data{115, 7}{1, 1}{1, 1}  = xpoints;
IC50data{115, 7}{1, 1}{1, 2}  = ypoints;
IC50data{115, 7}{3, 1}{1, 1}  = xpoints;
IC50data{115, 7}{3, 1}{1, 2}  = ypoints;
% clear dose response xpoints ypoints
%% plot IC50 curves
% plots mean IC50 curves for each cell line in three subplots corresponding
% to each drug. I plot the mean for each cell line and the mean of each
% individual fit, plus the mean IC50 estimated for each line
% I plot the BW control as there are no ancestor IC curves in triplicates
myDrugs = {'Streptozotocin' 'Triclabendazole' 'Sertraline HCl'};
myLines = {'L1' 'L2' 'L3' 'L4' 'L5'};
myColormap = [236 0 140]./255; 
% myColormap = [39 170 225]./255;
myColormap = [myColormap; myColormap; myColormap; myColormap; myColormap]; 
mSize=10; markerArray = {'d' 's' '^' '>' 'v'};
f=figure('Color','white');
for iDrug = 1:numel(myDrugs) %loop th
    curDrug = myDrugs{iDrug};
    curIC50 = IC50data(strcmp(IC50data(:,1),curDrug),:);
    IC50 = [];
    for iLine = 1:numel(myLines)
        curLine = myLines{iLine};
        myIC50 = curIC50(contains(curIC50(:,5),curDrug) & contains(curIC50(:,5),curLine),:);
        junk = myIC50{:,4}(1);
        junk = junk{:}(1);
        curDose = junk{:};
        curResponse = []; fitDose =[]; fitResponse=[];
        for i =1:3
            junk = myIC50{:,4}(i);
            junk = junk{:}(2);
            curResponse = [curResponse,junk{:}];
            junk = myIC50{:,7}(i);
            junk = junk{:}(1);
            fitDose =[fitDose; junk{:}];
            junk = myIC50{:,7}(i);
            junk = junk{:}(2);
            fitResponse =[fitResponse; junk{:}];
        end
        IC50 = [IC50;[mean(myIC50{1,3}(:,1)) mean(myIC50{1,3}(:,2))]];
        subplot(1,3,iDrug);hold on;
%         errorbar(log10(curDose),mean(curResponse,2),std(curResponse,[],2),markerArray{iLine},'MarkerFaceColor',myColormap(iLine,:),'MarkerEdgeColor',myColormap(iLine,:),'Color',myColormap(iLine,:),'MarkerSize',mSize);
        errorbar((curDose),mean(curResponse,2),std(curResponse,[],2),markerArray{iLine},'MarkerFaceColor',myColormap(iLine,:),'MarkerEdgeColor',myColormap(iLine,:),'Color',myColormap(iLine,:),'MarkerSize',mSize);
        plot((mean(fitDose)),mean(fitResponse),'Color',myColormap(iLine,:),'MarkerSize',mSize);
        %         errorbar(mean(IC50(:,1)),mean(IC50(:,2)), std(IC50(:,2)), std(IC50(:,2)),std(IC50(:,1)),std(IC50(:,1)),'o','MarkerFaceColor',myColormap(iLine,:),'MarkerEdgeColor',myColormap(iLine,:),'Color',myColormap(iLine,:))
%         plot(log10(mean(fitDose)),mean(fitResponse),'Color',myColormap(iLine,:),'MarkerSize',mSize);
        grid on; box on; ylim([0 1.2]); xlabel(['dose (uM)']);ylabel('normalized growth');
%         if iDrug==1
%             xlim([1 3])
%         elseif iDrug==2
%             xlim([0.5 2.5])
%         else
%             xlim([0 2])
%         end
    end
    myIC50 = curIC50(contains(curIC50(:,5),'BW wt'),:);
    junk = myIC50{1,4}(1);
    junk = junk{:}(1);
    curDose = junk{:};
    curResponse = []; fitDose =[]; fitResponse=[]; IC50bw=[];
    for i =1:3
        junk = myIC50{1,4}(i);
        junk = junk{:}(2);
        curResponse = [curResponse,junk{:}];
        junk = myIC50{1,7}(i);
        junk = junk{:}(1);
        fitDose =[fitDose; junk{:}];
        junk = myIC50{1,7}(i);
        junk = junk{:}(2);
        fitResponse =[fitResponse; junk{:}];
    end
    IC50bw =[IC50bw;[mean(myIC50{1,3}(:,1)) mean(myIC50{1,3}(:,2))]];
%     errorbar(log10(curDose),mean(curResponse,2),std(curResponse,[],2),'o','MarkerFaceColor','k','MarkerEdgeColor','k','Color','k','MarkerSize',mSize);
    errorbar(curDose,mean(curResponse,2),std(curResponse,[],2),'o','MarkerFaceColor','k','MarkerEdgeColor','k','Color','k','MarkerSize',mSize);
    %         errorbar(log10(mean(IC50(:,1))),mean(IC50(:,2)), std(IC50(:,2)), std(IC50(:,2)),std(IC50(:,1)),std(IC50(:,1)),'o','MarkerFaceColor','k','MarkerEdgeColor','k','Color','k')
%     plot(log10(mean(IC50bw(:,1))),mean(IC50bw(:,2)),'o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',12)
%     plot(log10(mean(fitDose)),mean(fitResponse),'Color','k');
%     plot(log10(mean(IC50(:,1))),mean(IC50(:,2)),'o','MarkerFaceColor',myColormap(iLine,:),'MarkerEdgeColor',myColormap(iLine,:),'MarkerSize',12)
   plot(mean(IC50bw(:,1)),mean(IC50bw(:,2)),'o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',12)
    plot(mean(fitDose),mean(fitResponse),'Color','k');
    plot(mean(IC50(:,1)),mean(IC50(:,2)),'o','MarkerFaceColor',myColormap(iLine,:),'MarkerEdgeColor',myColormap(iLine,:),'MarkerSize',12)
    set(gca, 'Xscale', 'log')
            if iDrug==1
            xlim([10^0.5 10^2.5])
        elseif iDrug==2
            xlim([10^0.5 10^2])
        else
            xlim([10^1 10^2.5])
        end

end
set(gcf,'Position',[1 595 1512 271])
% saveas(f,'ICcurves_fig5b_blue.svg')
% saveas(f,'ICcurves_fig5b_pink.svg')
%%