function montageOverlay(bv,fv,Alpha,dim,R,fig)
if ~exist('fig','var');fig =1;end
if ~exist('dim','var');dim =3;end
if ~exist('R','var');R =10;end


figure(fig)
clf;
ax1 = axes('Position',[0 0 1 1]);axis (ax1,'off');
ax2 = axes('Position',[0 0 1 1]);axis (ax2,'off');

imb = get.montage(bv,dim,R);
imagesc(ax1,imb);
colormap(ax1,'gray');
axis(ax1,'image');
if numel(size(fv))==3
    imf = get.montage(fv,dim,R);
else
    l = convert.dummy2label(fv);
    imf = get.montage(l,dim,R);
end
imagesc(ax2,imf,'AlphaData',(imf>0)*Alpha);
colormap(ax2,'jet');
axis(ax2,'off')
axis(ax2,'image')
linkaxes([ax1,ax2],'xy');


end