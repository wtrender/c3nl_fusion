function L = watershed(Y,r,ll)
% takes an intensity image and segments it using flood-fill 
if ~exist('r','var');r=1;end
if ~exist('ll','var');ll=0;end
L = ws(Y,r,1,ll);
end



% function L = blockWS(Y,dim,r,ll)
% tic
% g = ceil(dim./2^7);% subdivide large matrices to 2^7 parts
% N = numel(g);
% V = 1:max(g);
% partition = {};
% for s =1:numel(g)
%     partition{s} = [1:2^7:dim(s),dim(s)+1];
% end
% [t{N:-1:1}] = ndgrid(1:numel(V)) ;
% I = reshape(cat(N+1,t{:}),[],N) ;
% M = V(I);
% for s=1:N
% M(M(:,s)>g(s),:)=[];    
% end
% L = zeros(size(Y));
% m=1;
% for p=1:numel(M)
%     ind = cell(N,1);
%     for d =1:N
%         ind{d} = partition{d}(M(p,d)):partition{d}(M(p,d)+1)-1;
%     end
%     [L(ind{:}),m] = ws(Y(ind{:}),r,m);
%     fprintf('*')
%     toc
%     m
% end
% 
% end

function [L,m] = ws(Y,r,m,ll)
dim = size(Y);
if ~exist('m','var');m=1;end
%% the simple watershed algo 
if all(Y(:)==0 | Y(:)==1);% if binary convert to inverted distance matrix
    D = bwdist(~Y);
    D(~Y) = nan;
    Y = D;
    clear D;
end
P = padarray(Y,ones(numel(dim),1)*r,nan); % pad image to account for bounderies
Pdim =size(P);
d = c3nl.neighb(Pdim,r);
ix =find(~isnan(P));
[~,idx]=sort(-P(ix)); 
idx=ix(idx);
N=size(idx,1);
C=zeros(Pdim);
for ii=1:N
    ix = idx(ii)+d;
    c=C(ix);
    c=c(c>0);
    if isempty(c);C(idx(ii))=m;m=m+1;
    elseif ~any(diff(c));C(idx(ii))=c(1);
    else ;C(idx(ii))=mode(c);
    end
end
% plot.imOverlay(get.montage(Y,3),get.montage(C,3),0.4)
% stats = struct2array(regionprops(L,{'area'}))';
% ix = find(stats<=ll);
% d1 = c3nl.neighb(dim,r+1);% merge 1 volume class label with nearest connected labels
% tic
% for ii=1:numel(ix)
%    idx1=find(L==ix(ii));
%    idx2=[idx1+d1];
%    idx2 =idx2(~ismember(idx2(:),idx1));
%    idx = [idx1(:);idx2(:)];
%    if nnz((L(idx)>0))>2
%        [I,J,K]= ind2sub(dim,idx);
%        v = Y(idx);v(isnan(v))=0;
%        out = squareform(pdist([I',J',K']));
%        IDX  = ml.cluster.oKmeans([c3nl.scale(v),c3nl.scale(out,0,1,'col')/0.5],numel(unique(L(idx))));
%        L(idx1) = mode(L(idx(IDX==mode(IDX(1:numel(idx1))))));
%    end
% end
%C = convert.labelLUT(C);
stats = struct2array(regionprops(C,{'area'}))';
ix = find(stats<=ll);

for ii=1:numel(ix)
   idx1=find(C==ix(ii));
   idx2=[idx1'+d];
   idx2 =idx2(~ismember(idx2(:),idx1));
   [~,id]=min(abs(double(P(idx2(:)))-double(P(idx1(1)))));
   C(idx1) = C(idx2(id));
end
L = C;
for ii=1:numel(dim)
    switch ii
        case 1;L = L(1+r:end-r,:,:);
        case 2;L = L(:,1+r:end-r,:);
        case 3;L = L(:,:,1+r:end-r);
    end
end
L = atlas.reLabel(L);

end


