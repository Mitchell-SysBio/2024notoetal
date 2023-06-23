%% non-antibioitc cross-resistance analysis
% triclabendazole evolved lines were grown in the presence of lamotrigine
% or pyrimethamine with or without thymine to see if the folic acid
% synthesis pathway is also being targeted by this drugs (in addition to infb2)
% data from CL
% 04-26-2023 - MN
%% load data & process
load PYR_TCBZ_LAM_RoomTempIC50.mat
IC50_data = [];
for i = 1:length(IC50_total)
    curIC = IC50_total{:,i};
    IC50_data = cat(1,IC50_data,curIC);
end

IC50_data.Drug = {'Pyrimethamine' 'Pyrimethamine' 'Pyrimethamine' 'Pyrimethamine'...
    'Triclabendazole' 'Triclabendazole' 'Triclabendazole' 'Triclabendazole'...
    'Lamotrigine' 'Lamotrigine' 'Lamotrigine' 'Lamotrigine'}';
drugNames = unique(IC50_data.Drug);
%%
timePoint = 2; % 12 h
mSize = 8;
f=figure('Color','white');
for i = 1:numel(drugNames)
    curDrug = drugNames{i};
    curWT = [IC50_data.DoseResponse{contains(IC50_data.Drug,curDrug) & strcmp(IC50_data.FullName,'BW wt')}{timePoint}{:}];
    curEvo = [IC50_data.DoseResponse{contains(IC50_data.Drug,curDrug) & strcmp(IC50_data.FullName,'TCBZ Evo Line 2')}{timePoint}{:}];
    lineWT = [IC50_data.PredictedCurve{contains(IC50_data.Drug,curDrug) & strcmp(IC50_data.FullName,'BW wt')}{timePoint}{:}];
    lineEvo = [IC50_data.PredictedCurve{contains(IC50_data.Drug,curDrug) & strcmp(IC50_data.FullName,'TCBZ Evo Line 2')}{timePoint}{:}];
    stdWT = [IC50_data.Stds{contains(IC50_data.Drug,curDrug) & strcmp(IC50_data.FullName,'BW wt'),timePoint}{:}];
    stdEvo = [IC50_data.Stds{contains(IC50_data.Drug,curDrug) & strcmp(IC50_data.FullName,'TCBZ Evo Line 2'),timePoint}{:}];
    subplot(1,3,i); hold on;
    errorbar(curWT(:,1),curWT(:,2),stdWT,'o','MarkerFaceColor','k','MarkerEdgeColor','k','Color','k','MarkerSize',mSize);
    plot(lineWT(1:5000),lineWT(5001:end),'k');
    errorbar(curEvo(:,1),curEvo(:,2),stdEvo,'o','MarkerFaceColor','m','MarkerEdgeColor','m','Color','m','MarkerSize',mSize);
    plot(lineEvo(1:5000),lineEvo(5001:end),'m');
    xlim([1 550]); ylim([0 1.2]); title(curDrug);
    if i==3
         xlim([0.3 550]); ylim([0 1.2]); title(curDrug);
    end
    grid on; box on; 
    set(gca, 'Xscale', 'log'); 
end
set(gcf,'position',[-21        1265        1057         223]);
% 
% saveas(f,'nabxICcurves.fig');
% saveas(f,'nabxICcurves.svg');
