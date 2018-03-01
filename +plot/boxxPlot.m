function boxxPlot(X,gp,cmap,order,fuzzy,w,fig,ax,ms,xl)
% X is a 1d long vector that corropsonds to hierarchy in gp


gp = categorical(gp);
gn = categories(gp);
if ~exist('gap','var');gap = 0.2;end
if ~exist('ms','var');ms = 2;end
if ~exist('cmap','var');cmap = jet(d(2)*d(3));end
if ~exist('fuzzy','var');fuzzy = 0;end
if ~exist('lim','var');lim = [min(X(:))*1.005,max(X(:))*1.005];end
if ~exist('fig','var');figure(1);else ;figure(fig);end
if ~exist('w','var');w = 0.35;end
if ~exist('xl','var');xl = 1;end
if ~exist('ax','var')|isempty(ax);ax = gca;end


% start by getting the quartiles 
Q = zeros(numel(order),10);% first,median,third,IQR,lower inner fence, lower outer fence ,upper inner fence ,upper outer fence,lower whiskers,upper whiskers 
for ii=1:numel(order)
    tmp = X(gp==order{ii});
    Q(ii,1:3) = quantile(tmp,[0.25 0.50 0.75]);
    Q(ii,4) = range(Q(ii,1:3));
    Q(ii,5:8) = [Q(ii,1)-Q(ii,4)*1.5,Q(ii,1)-Q(ii,4)*3,Q(ii,3)+Q(ii,4)*1.5,Q(ii,3)+Q(ii,4)*3];
    Q(ii,9:10) = [min(tmp(tmp>=Q(ii,5))),max(tmp(tmp<=Q(ii,7)))];
end

hold all;

for c=1:numel(order)  
    tmp = X(gp==order{c});
    P = plot.bezierrect(c-w/2,Q(c,1),w,Q(c,4),0.05); % get curved box patch
    if fuzzy
        plot(c*(tmp.^0)+c3nl.scale(randn(size(tmp)),-0.25,0.25),tmp,'Marker','o','LineStyle','none','MarkerFaceColor',cmap(c,:),'MarkerEdgeColor','none','MarkerSize',ms);  
    end
    patch(P(:,1),P(:,2),cmap(c,:).^.5,'FaceAlpha',0.5,'parent',ax,'edgecolor','none');
    line(P(:,1),P(:,2),'linewidth',1,'color','k','parent',ax)
    line([c-w/2,c+w/2],[Q(c,2),Q(c,2)],'linewidth',1,'color','k','parent',ax);
    line([c,c],[Q(c,3),Q(c,10)],'linewidth',1,'color','k','parent',ax);
    line([c-w/4,c+w/4],[Q(c,10),Q(c,10)],'linewidth',1,'color','k','parent',ax);
    line([c,c],[Q(c,1),Q(c,9)],'linewidth',1,'color','k','parent',ax);
    line([c-w/4,c+w/4],[Q(c,9),Q(c,9)],'linewidth',1,'color','k','parent',ax);

end
axis([0.5 numel(order)+0.5 ylim])
ax.XTick = 1:numel(order);
if xl
    ax.XTickLabel = order;
    ax.XTickLabelRotation=-45;
else 
    ax.XTickLabel = [];
end


end