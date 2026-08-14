function percentiles=getPercentiles(data,percentileSet)
percentiles = [];
if(~isempty(data))
    assert(min(percentileSet)>= 0 && max(percentileSet) <= 1);
    
    n=percentileSet*length(data);
    n=round(n);
    dataOrder=sortrows(data);
    
    n(n<=0)=1;
    n(n>length(dataOrder))=length(dataOrder);

    percentiles = dataOrder(n);
end
