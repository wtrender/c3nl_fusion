function covWithWeights(X,W,cmap)
disp('')
d = size(X);
[x,y,h,w]=gen.fancyGrid(d(2)*[1,1],0.05,0.05,'weighted',[0.7,0.3]);
figure();
clf
mx = 0;mn=0;
for ii=1:d(2)
    r= ii+d(2);
    axh.(['P' num2str(ii)]) = axes('position',[x(ii),y(ii),h(ii),w(ii)]);
    im = cov(X{ii});
    imagesc(im,'parent',axh.(['P' num2str(ii)]));
if ~isempty(W)
    axh.(['P' num2str(r)]) = axes('position',[x(r),y(r),h(r),w(r)]);
    ww = W{ii};
    if mx<max(ww);mx=max(ww);end
    if mn>min(ww);mn=min(ww);end
for jj=1:size(im,1)
    [xx,yy]=plot.square(jj,ww(jj)/2,0.5,ww(jj));
    patch(xx,yy,cmap(ii,:),'FaceAlpha',0.75,'parent',axh.(['P' num2str(r)]),'Edgecolor','k','linewidth',0.5);
    
end
axh.(['P' num2str(r)]).XTick = 1:size(im,1);
axh.(['P' num2str(r)]).XTickLabel = '';
end
axh.(['P' num2str(ii)]).YTick = 1:size(im,1);
axh.(['P' num2str(ii)]).XTickLabel = '';
end
if ~isempty(W)
for ii=1:d(2)
   axis(axh.(['P' num2str(ii+d(2))]),[0.5,numel(ww)+0.5,min(mn+sign(mn)*mx*0.1,0),mx*1.1])
end
end

end