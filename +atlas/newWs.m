function [L,stat,T] = newWs(Y,r,ex,ll,ul)
if ~exist('r','var')||isempty(r);r=1;end
if ~exist('ex','var')||isempty(ex);ex=0;end
if all(Y(:)==0 | Y(:)==1);% if binary convert to inverted distance matrix
    D = bwdist(~Y);
    D(~Y) = nan;
    Y = D;
    clear D;
end
CC = bwconncomp(~isnan(Y),conndef(ndims(Y),'maximal'));
ix = cellfun(@numel, CC.PixelIdxList)<ex;
% exlude small volumes  
for ii=1:numel(ix)
    if ix(ii)
        Y(CC.PixelIdxList{ii})=nan;
    end
end
tic
CC = bwconncomp(~isnan(Y));
LL=zeros(size(Y));TT=zeros(size(Y));
for jj=1:CC.NumObjects
y = -inf(size(Y));
y(CC.PixelIdxList{jj})=Y(CC.PixelIdxList{jj});
Pdim =size(y);
d = c3nl.neighb(Pdim,r); % get connected components linear index
[~,idx] = sort(-y(:)); % force minima
N=size(idx,1);
L = zeros(Pdim);
T = zeros(Pdim);
m=1;
mxv = numel(y);
for ii=1:N
    idx1=idx(ii);
    if y(idx1) == -inf;break;end
    idx2=idx1+d;
    c = L(idx2(idx2>0&idx2<mxv));% get voxel neighborhood
    l  = sort(c(:)); 
    l(~diff(l)) = [];
    l(l==0 ) = [];
    t = find([isempty(l),numel(l)==1,numel(l)>1],1,'first');
    switch t
        case 1;L(idx(ii))=m;m=m+1;T(idx(ii))=0;% local minima with no surrounding labels
        case 2;L(idx(ii))=l(1);T(idx(ii))=1;% local minima with no surrounding labels
        case 3        
        L(idx(ii))=mode(c(c>0));% Majority vote wins and makes smooth parcels
        T(idx(ii))=3;
    end
    if ~mod(ii,1000)
        if mod(ii,50000);fprintf('*');else;fprintf('\n');end
    end
end
LL = atlas.reLabel(atlas.merge(double(LL),L*100));
TT = atlas.merge(TT,T);
end
L = LL;
T = TT;
%% merge below lower limit parcels with nearest cluster
if exist('ll','var')&&~isempty(ll)
d2 = c3nl.neighb(Pdim,1,'city');
ix=1;
while ~isempty(ix)
    L = atlas.reLabel(L);
    stat = regionprops3(L,Y,{'Volume','MeanIntensity','Centroid'});
    stat.id = (1:height(stat))';
    stat = sortrows(stat,{'Volume'});
    ix = find(stat.Volume<=ll);
    for ii=1:numel(ix)
       idx1=find(L==stat.id(ix(ii)));
       idx2=idx1'+d2;
       idx2 = idx2(idx2>0&idx2<mxv);
       idx2 = idx2(~ismember(idx2(:),idx1));
       [~,lut]=convert.label2Dummy(L(idx2(:)));
       if ~isempty(lut)
        [~,idx3]=min(dist.euc(stat{ix(ii),:},stat{lut,:}));
        L(idx1) = lut(idx3);
       else 
           L(idx1) = 0;
       end
    end
end
L = double(atlas.reLabel(L));
end

% split large parcels to upper limit parcels
if exist('ul','var')&&~isempty(ul)
stat = struct2array(regionprops(L,{'area'}))';
ix = find(stat>ul);
m = max(L(:));
for ii=1:numel(ix)
   idx1=find(L==ix(ii));
   [s{1:numel(dim)}] = ind2sub(Pdim,idx1(:));
   k = ceil(stat(ix(ii))/ul);
   id= kmeans([c3nl.scale(cell2mat(s),0,1,'col'),Y(idx1)],k);
   L(L==ix(ii)) = m+id;
   m=max(L(:));
end
L = atlas.reLabel(L);
end
stat = regionprops3(L,Y,{'Volume','MeanIntensity'});

L(Y==-inf)=nan;
toc



end