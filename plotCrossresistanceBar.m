load('EvoLinesIC50Results.mat');
foldChangeTable2 = foldChangeTable;
junk = find(contains(foldChangeTable2.strainList,'Sertra'));
junk(21:25,:)=[];
foldChangeTable2(junk,:)=[];
load('EvoLinesIC50Results2.mat');
foldChangeTable = cat(1,foldChangeTable,foldChangeTable2);
foldChangeTable.drugList(startsWith(foldChangeTable.drugList,'Dox')) = {'Doxycycline'};
myDrugs = {'Streptozotocin' 'Triclabendazole' 'Sertraline HCl'};
abDrugs = {'Doxycycline' 'Erythromycin' 'Ciprofloxacin' 'Carbenicillin' 'Sulfamonomethoxine' 'Furazolidone' };
myLines = {'L1' 'L2' 'L3' 'L4' 'L5'};
mechColor = [186,85,211; 125,25,207; 220,20,60; 255,215,0; 0,128,0;240,128,128];
mechColor = mechColor./255;
% ratiosAB = [4.8/0.51 54.17/14.7 0.16/0.04 37.6/17 0.28/0.13 18.73/2.7];


%% plot bar
f=figure('Color','white');
for iDrug = 1:numel(myDrugs)
    curDrug = myDrugs{iDrug};
    curDrugData = foldChangeTable(contains(foldChangeTable.strainList,curDrug),:);
    for iAb = 1:numel(abDrugs)
        curAb = abDrugs{iAb};
        curFC = curDrugData{strcmp(curDrugData.drugList,curAb),3};
%         curFC(isnan(curFC)) = ratiosAB(iAb);
        curFC = log2(curFC);
        subplot(1,3,iDrug); hold on;
        %         plot(linspace(iAb,iAb,numel(curFC))',curFC,'o','MarkerEdgeColor',mechColor(iAb,:),'MarkerFaceColor',mechColor(iAb,:))
        b = bar(iAb,median(curFC),'FaceColor','flat');
        b.CData = mechColor(iAb,:);
        plot(linspace(iAb,iAb,numel(curFC))',curFC,'o','MarkerEdgeColor','k','MarkerFaceColor','k')

        %         legend(abDrugs)
    end
    grid on; box on; xlim([0 7]); ylim([-log2(64) log2(64)]) %ylim([0 7])
    xticks([1:6]);xticklabels(abDrugs); xtickangle(45)
    ylabel('log2(IC50 fold change)');
    yline(0,'r')

end
set(gcf,'Position',[1 541 1512 325])
