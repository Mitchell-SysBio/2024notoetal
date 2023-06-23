function [pwNames, pvalFDR, pwGenes, pwGenesUniv pwGenesSet]= checkPwEnrichment (geneSet,geneUniverse,pwAnn,pwSizeRange)

if (nargin==3)
    pwSizeRange = [10 1000];
end

pwGenes = {};pwGenesUniv = {};pwGenesSet = {}; % added 02-08-2023
pwSizes = cellfun(@numel,{pwAnn.genes})';
p=[];
k=1;
for iPw = 1:length(pwAnn)
    if (pwSizes(iPw) > pwSizeRange(1,1) & pwSizes(iPw) < pwSizeRange(1,2))
        curPw = pwAnn(iPw).genes;
        if (length(intersect(curPw,geneSet))>1)
            if (length(intersect(geneUniverse,geneSet)) < length(geneSet))
                disp 'error: gene set not contained in gene universe'
            end
            g1 = geneSet;
            g2 = intersect(geneUniverse,curPw);
            N_nono = length(setdiff(geneUniverse,cat(1,g1,g2)));
            N_yesyes = length(intersect(g1,g2));
            N_yesno = length(setdiff(g1,g2));
            N_noyes = length(setdiff(g2,g1));
            T = [N_nono, N_noyes; N_yesno, N_yesyes];
            [d(k),p(k)] = fishertest(T); % d = decision, p = pvalue
            testedPW(k)=iPw;
            pwGenes{k} = intersect(g1,g2);  % added 02-08-2023
            pwGenesUniv{k} = g2;  % added 02-08-2023
            pwGenesSet{k} = curPw;
            k=k+1;
        end
    end
end
if (length(p)==0)
    pvalFDR = nan;
    pwNames = 'no tests were made';
    genes = 'no tests were made'; % added 02-08-2023
else
    pvalFDR = mafdr(p,'bhfdr',1)';
pwNames = {pwAnn(testedPW).names}';
% genes = intersect(g1,g2); % added 02-08-2023
end


end

