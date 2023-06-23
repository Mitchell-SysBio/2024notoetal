function [pval] = calcEmpiricalPVal(data,observation,nShuffles)

if(nargin==2)
    nShuffles = 100000;
end

inx = randi(length(data),nShuffles,length(observation)); % get random indexes (sized by nShuffles and observation size)
shuffledData = nan(nShuffles,length(observation)); % place holder for shuffled data
shuffledData(:) = data(inx(:)); % populate shuffled data

pval = sum(nanmean(shuffledData,2)<=nanmean(nanmean(observation,2)))/nShuffles;

end