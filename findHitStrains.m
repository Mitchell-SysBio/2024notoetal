function [tf pvalFDR] = findHitStrains(LFC,inxDrugs,tfClusterDrugs,pvalFDRCutoff)

if(nargin==3), pvalFDRCutoff = 0.05; end

% before anything, impute missing values by k nearest rows (strains)
data = LFC';
imputedData = knnimpute(data,5);
LFC2 = imputedData';

ia = inxDrugs(tfClusterDrugs); ib = inxDrugs(~tfClusterDrugs); % find the indexes of both drug groups
% itergate over genes and ttest each
for i=1:size(LFC2,1)
    [h pval(i)] = ttest2(LFC2(i,ismember(inxDrugs,ia)),LFC2(i,ismember(inxDrugs,ib)));
end

[pvalFDR] = mafdr(pval,'bhfdr',1);
tf = false(size(pvalFDR));
tf(pvalFDR<pvalFDRCutoff) = true;

end