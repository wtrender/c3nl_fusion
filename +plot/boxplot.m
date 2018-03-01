function boxplot(X,labels,mode,cmap,fc,fs,fig,lim,yx,ylabel)
disp('');
d = size(X);
if numel(d)<3;d(3) =1;end
if ~exist('gap','var');gap = 0.2;end
if ~exist('cmap','var');cmap = jet(d(2)*d(3));end
if ~exist('fuzzy','var');fuzzy = 0;end
if ~exist('lim','var');lim = [min(X(:))*1.005,max(X(:))*1.005];end
if ~exist('yx','var');yx = 0;end
if ~exist('fig','var');figure(1);else ;figure(fig);end
if ~exist('fc','var');fc=0.02;end


% start by getting the quartiles 
Q = zeros(prod(d(2:3)),10);% first,median,third,IQR,lower inner fence, lower outer fence ,upper inner fence ,upper outer fence,lower whiskers,upper whiskers 
for ii=1:prod(d(2:3))
    tmp = X(:,ii);
    Q(ii,1:3) = quantile(tmp,[0.25 0.50 0.75]);
    Q(ii,4) = range(Q(ii,1:3));
    Q(ii,5:8) = [Q(ii,1)-Q(ii,4)*1.5,Q(ii,1)-Q(ii,4)*3,Q(ii,3)+Q(ii,4)*1.5,Q(ii,3)+Q(ii,4)*3];
    Q(ii,9:10) = [min(tmp(tmp>=Q(ii,5))),max(tmp(tmp<=Q(ii,7)))];
end

[y,x,h,w] = get.fancyGrid(ones(d(3),1),0.01,0.2,'stretch');


clf;
ax = [];
for ii=1:d(3)
    ax.(['p' num2str(ii)]) = axes('position',[x(ii),y(ii),w,h]);
    grid on
    
end

hold all;
c = 1;

for ii=1:d(3)   
    w = .5;
    mny = 0;mxy =0;
    for jj=1:d(2)
        xt = X(:,c);
        if max(xt(:))>lim(2);lim(2)=max(xt(:))+eps;end
        if min(xt(:))<lim(1);lim(1)=min(xt(:))-eps;end
        P = plot.bezierrect(jj-w/2,Q(c,1),w,Q(c,4),0.05); % get curved box patch
        yt = X(:,c);
        patch(P(:,1),P(:,2),cmap(jj,:).^.5,'FaceAlpha',0.5,'parent',ax.(['p' num2str(ii)]),'edgecolor','none');
        line(P(:,1),P(:,2),'linewidth',0.5,'color','k','parent',ax.(['p' num2str(ii)]))
        line([jj-w/2,jj+w/2],[Q(c,2),Q(c,2)],'linewidth',0.5,'color','k','parent',ax.(['p' num2str(ii)]));
        line([jj,jj],[Q(c,3),Q(c,10)],'linewidth',0.5,'color','k','parent',ax.(['p' num2str(ii)]));
        line([jj-w/4,jj+w/4],[Q(c,10),Q(c,10)],'linewidth',0.5,'color','k','parent',ax.(['p' num2str(ii)]));
        line([jj,jj],[Q(c,1),Q(c,9)],'linewidth',0.5,'color','k','parent',ax.(['p' num2str(ii)]));
        line([jj-w/4,jj+w/4],[Q(c,9),Q(c,9)],'linewidth',0.5,'color','k','parent',ax.(['p' num2str(ii)]));
        c = c+1;
    end
end
for jj=1:d(3)
    %axis(ax.(['p' num2str(ii)]),[0,1,0,maxdens*1.1]);
    grid on
    axis(ax.(['p' num2str(jj)]),[0.25,d(2)+.75,lim]);
    if jj>1||~yx       
        ax.(['p' num2str(jj)]).YTickLabel ='';        
    end
    ax.(['p' num2str(jj)]).XTick = 0:d(3);
    ax.(['p' num2str(jj)]).XTickLabel ='';
    ax.(['p' num2str(jj)]).Color = ones(1,3)*0.90;  
    if exist('labels','var')&& ~isempty(labels)
    title(ax.(['p' num2str(jj)]),labels{2}{jj})
    end
    if jj==1&yx 
        ax.p1.YLabel.String = ylabel;
    end
%     x = xlim;
%     text(ax.(['p' num2str(jj)]),range(x)/2,-1,labels{1,1}{jj});
    
    %ax.(['p' num2str(ii)]).YLabel.Position = [102.5,ax.(['p' num2str(ii)]).YLabel.Position(2:3)];
end


end


