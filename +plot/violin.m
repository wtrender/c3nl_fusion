function violin(X,labels,cmap,fig,lim,yx,ylabel,fuzzy)
disp('');
d = size(X);
if numel(d)<3;d(3) =1;end
if ~exist('gap','var');gap = 0.2;end
if ~exist('cmap','var');cmap = jet(d(2)*d(3));end
if ~exist('fuzzy','var');fuzzy = 0;end
if ~exist('lim','var');lim = [min(X(:))*1.005,max(X(:))*1.005];end
if ~exist('yx','var');yx = 0;end


% start by getting the quartiles 
Q = zeros(prod(d(2:3)),10);% first,median,third,IQR,lower inner fence, lower outer fence ,upper inner fence ,upper outer fence,lower whiskers,upper whiskers 
for ii=1:prod(d(2:3))
    tmp = X(:,ii);
    Q(ii,1:3) = quantile(tmp,[0.25 0.50 0.75]);
    Q(ii,4) = range(Q(ii,1:3));
    Q(ii,5:8) = [Q(ii,1)-Q(ii,4)*1.5,Q(ii,1)-Q(ii,4)*3,Q(ii,3)+Q(ii,4)*1.5,Q(ii,3)+Q(ii,4)*3];
    Q(ii,9:10) = [min(tmp(tmp>=Q(ii,5))),max(tmp(tmp<=Q(ii,7)))];
end

[y,x,h,w] = get.fancyGrid(ones(d(2),1),0.01,0.2,'stretch');
%x = x*100;w=w*100;
if ~exist('fig','var');figure(1);else ;figure(fig);end
if ~exist('fc','var');fc=0.02;end

clf;
ax = [];
for ii=1:d(2)
    ax.(['p' num2str(ii)]) = axes('position',[x(ii),y(ii),w,h]);
    grid on
    
end

hold all;
c = 0;
for ii=1:d(3)   
    for jj=d(2):-1:1
        c = c+1;
        xt = X(:,jj,ii);
        if max(xt(:))>lim(2);lim(2)=max(xt(:))+min(xt(:))*0.1;end
        if min(xt(:))<lim(1);lim(1)=min(xt(:))-min(xt(:))*0.1;end
        [f,xi] = ksdensity(xt,'Support',lim);
        f = (f./(max(f)*2.5));
        f(f<1e-3)=0;
        %[a,b]=hist(xt,20);
        %plot3(xi,ones(numel(xi),1)*jj,f,'color',cmap(c,:));
        patch([f,fliplr(-f)]+ii,[xi,fliplr(xi)],cmap(ii,:),'FaceAlpha',0.5,'parent',ax.(['p' num2str(jj)]),'Edgecolor',max(cmap(ii,:)-0.5,0));
        if fuzzy
            for k=1:numel(xi)-1
                tmp = xt(xt>xi(k)&xt<=xi(k+1));
                if ~isempty(tmp)
                    xx = linspace(ii-f(k),ii+f(k),numel(tmp)+2);
                    for p=1:numel(tmp)
                       hold(ax.(['p' num2str(jj)]), 'on')
                       plot(xx(p+1),tmp(p),'marker','o','parent',ax.(['p' num2str(jj)]),'markerfacecolor',max(cmap(ii,:)-0.3,0),'markersize',fuzzy,'markeredgecolor','none');
                       hold(ax.(['p' num2str(jj)]), 'off')
                    end
                end
            end
        end
        if min(ylim)<0
            line(xlim,[0,0],'LineStyle','--','Color','k','parent',ax.(['p' num2str(jj)]));
        end
        %patch(-f+ii*0.2,xi,cmap(ii,:),'FaceAlpha',0.5,'parent',ax.(['p' num2str(jj)]));

    end
end

for jj=1:d(2)
    %axis(ax.(['p' num2str(ii)]),[0,1,0,maxdens*1.1]);
    grid on
    axis(ax.(['p' num2str(jj)]),[0,d(3)+1,lim]);
    if jj>1||~yx       
        ax.(['p' num2str(jj)]).YTickLabel ='';        
    end
    ax.(['p' num2str(jj)]).XTick = 0:d(3);
    ax.(['p' num2str(jj)]).XTickLabel ='';
    ax.(['p' num2str(jj)]).Color = ones(1,3)*0.90;  
    if exist('labels','var')&& ~isempty(labels)
    title(ax.(['p' num2str(jj)]),labels{1}{jj})
    end
    if jj==1&yx 
        ax.p1.YLabel.String = ylabel;
    end
%     x = xlim;
%     text(ax.(['p' num2str(jj)]),range(x)/2,-1,labels{1,1}{jj});
    
    %ax.(['p' num2str(ii)]).YLabel.Position = [102.5,ax.(['p' num2str(ii)]).YLabel.Position(2:3)];
end
if yx
lg = legend(ax.p1,labels{1,2},'Location','northeast','Orientation','horizontal');
lg.Position = [0.2,0.05,0.3,0.05];
lg.Box = 'off';
end

    

end