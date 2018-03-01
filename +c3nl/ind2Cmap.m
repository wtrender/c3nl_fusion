function c = ind2Cmap(map,v,n)
% map = a matrix of at least two colors 
% v a vector of data values
if size(map,1)==1
    map = [0.5,0.5,0.5;map];
end
n=1000;
cmap = get.cmap(map,n);
ix = floor(1+ (v(:) - min(v(:))) / (range(v(:))) * (n-1));
if any(isnan(ix));
    c = zeros(numel(ix),3);
    c(isnan(ix),:) = repmat(map(end,:),nnz(isnan(ix)),1);
    c(~isnan(ix)) = cmap(~isnan(ix),:);
else
    c = cmap(ix,:);
end
end