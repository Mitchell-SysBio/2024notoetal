function [pwAnn] = extractPwAnnFromGMT(strFilename,minPwSize,maxPwSize)
if (nargin==1)
    minPwSize = 1;
    maxPwSize = 2000;
    disp 'max patwhway size will be set to 2000 by default'
end
strainMetadata = readcell('metadata/10_2022_strainMetadata.csv');
strainNames = strainMetadata(2:end,1);
strainEntrez = cell2mat(strainMetadata(2:end,20));

pwRaw = readcell(strFilename,'Delimiter',"\n");
pwAnn = struct;
k=1;
n=[];
for i=1:size(pwRaw,1)
    junk = strsplit(pwRaw{i},'\t');
    if (size(junk,2)>(minPwSize+1) & size(junk,2)<maxPwSize)
        pwAnn(k,1).names = junk{1};
        entrezList = cell2mat(cellfun(@str2num,junk(2:end),'UniformOutput',false));
        inxStrains = find(ismember(strainEntrez,entrezList));
         pwAnn(k).genes = strainNames(inxStrains);

        k=k+1;
    end
             n(i) = length(junk)-1; 
end
% figure;hold on; bar(sort(n));yline(minPwSize);yline(maxPwSize);title([strFilename ' - saved ' num2str(k-1) ' out of ' num2str(i) ' pathways'])
end