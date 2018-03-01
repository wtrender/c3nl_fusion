function labels(T,pxw,pxh,C,cmap)

figure('units','pixels','position',[0 0 pxw pxh],'Visible','on','resize','off');
clf
n = height(T);
k = ceil(n/C);
T = sortrows(T,'cid','descend');
[x,y] = get.fancyGrid(ones(1,k),0.04,0.03,'same');
[yc,xc,h,w] = get.fancyGrid(ones(C,1),0.04,0.03,'same');
labels = fliplr(strrep(T.name,'_',' '));
num = fliplr(T.number);
c=0;
for jj=C:-1:1
    ax = axes('position',[xc(jj) yc(jj) w h]);
    for ii=1:k
        if c<n
        c=c+1;
        text(x(ii),y(ii),{char(9679)},'FontUnits','normalized','fontsize',0.1,'HorizontalAlignment','center','parent',ax,'color',cmap(T.cid(c),:))
        text(x(ii),y(ii),sprintf('%i',num(c)),'FontUnits','normalized','fontsize',0.03,'HorizontalAlignment','center','FontWeight','bold','color','w','parent',ax);
        text(x(ii)+0.1,y(ii),labels{c},'HorizontalAlignment','left','FontUnits','normalized','fontsize',0.03,'FontWeight','bold','parent',ax);
        end
    end
    axis off
end
end
