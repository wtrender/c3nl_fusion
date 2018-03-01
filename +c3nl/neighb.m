function [D,W,w] = neighb(od,r,dist)

% get back linear shift for neighbourhood
% od = original matrix size
% pd = paddedd matrix size
% r = neighbourhood radius
if ~exist('dist','var');dist='euc';end
D = [];
[D{1:length(od)}]=ndgrid(-r:r);
d = cellfun(@(x) x+1+r,D,'un',0);
d=sub2ind(od,d{:});
d=d-d((numel(d)+1)/2);
switch dist
    case 'euc';w = (bwdist(d==0,'euclidean'));
    case 'city';w = (bwdist(d==0,'cityblock'));
    case 'chess';w = (bwdist(d==0,'chessboard'));
    case 'qeuc';w = (bwdist(d==0,'quasi-euclidean'));
end

D = d(w<=r*sqrt(2));
W = double(1./w(w<=r*sqrt(2)));W(D==0)=[];
D(D==0)=[];
w = (w>=0)./numel(w);



end