function histfit(X,labels,lim,cmap,fc,fig)
disp('');
[m,n]= size(X);
xx = X(:);
if ~exist('cmap','var');cmap = jet(m*n);end
% 
% start by getting the quartiles 
Q = zeros(m*n,10);% first,median,third,IQR,lower inner fence, lower outer fence ,upper inner fence ,upper outer fence,lower whiskers,upper whiskers 
pmap = rgb2hsv(cmap);
pmap(:,2) = 0.55;
pmap(:,3) = 0.75;
pmap = hsv2rgb(pmap);
for ii=1:m*n
    tmp = xx{ii}(:);
    if ~isempty(tmp)
        Q(ii,1:3) = quantile(tmp,[0.25 0.50 0.75]);
        Q(ii,4) = range(Q(ii,1:3));
        Q(ii,5:8) = [Q(ii,1)-Q(ii,4)*1.5,Q(ii,1)-Q(ii,4)*3,Q(ii,3)+Q(ii,4)*1.5,Q(ii,3)+Q(ii,4)*3];
        Q(ii,9:10) = [min(tmp(tmp>=Q(ii,5))),max(tmp(tmp<=Q(ii,7)))];
    end
end

[x,y,w,h] = get.fancyGrid(ones(m,1),0.04,0.1,'stretch');
%x = x*100;w=w*100;
if ~exist('fig','var');figure(1);else ;figure(fig);end
if ~exist('fc','var');fc=0.02;end

clf;
ax = [];
for ii=1:m
    ax.(['p' num2str(ii)]) = axes('position',[x(ii),y(ii)+0.1,w,h]);
    grid on
    
end

hold all;
c = 0;
maxdens = 0;
label = labels{1,2};
for ii=1:n   
    for jj=1:m
        xt = X{jj,ii}(:);
        if ~isempty(xt)
            [f,xi] = ksdensity(xt,'Support',lim);
            if maxdens<max(f);maxdens=max(f);end
            patch(xi,f,cmap(ii,:),'FaceAlpha',0.5,'parent',ax.(['p' num2str(jj)]));
        end
    end
end

for ii=1:m
    axis(ax.(['p' num2str(ii)]),[lim,ylim(ax.(['p' num2str(ii)]))]);
    grid on
    ax.(['p' num2str(ii)]).YTickLabel ='';
    ylabel(ax.(['p' num2str(ii)]),labels{1,1}{ii})
    %ax.(['p' num2str(ii)]).YLabel.Position = [102.5,ax.(['p' num2str(ii)]).YLabel.Position(2:3)];
    if ii>1;ax.(['p' num2str(ii)]).XTickLabel ='';end
end
% legend(ax.(['p' num2str(ii)]),labels{1,2});
% xlabel(ax.p1,'Accuracy')


end