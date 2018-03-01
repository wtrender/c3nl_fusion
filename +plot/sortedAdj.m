function sortedAdj(w,T,ix,cmap,clim,fig)
disp('')
n = height(T);
adj = zeros(n);
[~,b]=boundary(paretotails(w,0.025,0.975));
W = max(min(w,b(2)),b(1));
adj(triu(true(n),1)) = W;
adj = symAdj(adj,'upper');
if ~exist('fig','var');figure(1);else;figure(fig);end
if ~exist('ix','var');ix = 1:n;end
clf
ax = axes('Position',[0.2,.2,.6,.6]);
imagesc(adj(ix,ix));

W = sum(dummyvar(T.cid(ix)));
C = cumsum(W);

for ii=1:numel(C)
[x,y] = plot.square(C(ii)-W(ii)/2+.5,-.8,W(ii),2);
patch(x,y,cmap(ii,:),'parent',ax,'edgecolor','none');%,'edgecolor','none'
patch(y,x,cmap(ii,:),'parent',ax,'edgecolor','none');
end
axis(ax,[-1.5,C(end),-1.5,C(end)]);
daspect(ax,[1,1,1])
axis(ax,'off')
colormap(ax,parula)
if exist('clim','var');ax.CLim = clim;end
end