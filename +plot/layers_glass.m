function [img,mask,cmap] = layers_glass(vol,G,or,method,thr,T,cmap,bkg,fig)

load +atlas/cortex
dpi = 96;
ao = layers.(or).ao.img;
ao_alpha = layers.(or).ao_alpha.img;
mask = layers.(or).Alpha;
lim = layers.(or).ao.lim;
az = layers.(or).az;
el = layers.(or).el;
if ~exist('thr','var');thr=2;end
if ~exist('T','var');T=[];end
if ~exist('cmap','var');cmap=hot;end
if ~exist('bkg','var');bkg=[246,246,256]./256;end

switch method
    case 'mip'
        im = plot.roi2imag(vol,struct('G',G,'type','mip','az',az,'el',el,'cmap',cmap,'marker','.','thr',thr,'lim',lim,'dpi',dpi));
    case 'pc'
        im = plot.roi2imag(vol,struct('G',G,'type','pc','az',az,'el',el,'cmap',cmap,'marker','.','thr',thr,'lim',lim,'dpi',dpi));
    case 'roi'
        im = plot.roi2imag(vol,struct('G',G,'type','roi','az',az,'el',el,'cmap',cmap,'marker','.','thr',thr,'lim',lim,'dpi',dpi,'Smooth',0.6));
    case 'nodes'
        im = plot.roi2imag(vol,struct('G',G,'type','nodes','az',az,'el',el,'cmap',cmap,'T',T,'lim',lim,'dpi',dpi,'adj',eye(height(T))));
    case 'stack'
        im = plot.roi2imag(vol,struct('G',G,'type','stack','az',az,'el',el,'cmap',cmap,'T',T,'lim',lim,'dpi',dpi,'bkg',[.5,.5,.5]));    
end

d = size(ao);
if ~exist('fig','var');fig = 1;end
fg = figure(fig);clf;
set(fg,'units','pixels','position',[0 0 d(2) d(1)],'Visible','on','resize','off');
Alpha = 0.8;
fc = 'w';

ax0 = axes('position' ,[0.02 0.05 0.9 0.9]);
imagesc(mask,'AlphaData',mask*0.25,'parent',ax0);
axis(ax0, 'off');axis(ax0,'image');
Alpha = 0.5;
%fc = 'k';
colormap(ax0,[0.5,0.5,0.7])

ax1 = axes('position' ,[0.02 0.05 0.9 0.9]);
imagesc(im.img,'AlphaData',im.alpha.*mask,'parent',ax1);
axis(ax1, 'off');axis(ax1,'image');
ax2 = axes('position' ,[0.02 0.05 0.9 0.9]);
imagesc(ao,'AlphaData',c3nl_scale(double(ao_alpha).*(mask))*Alpha,'parent',ax2);
axis(ax2, 'off');axis(ax2,'image')
colormap(ax2,bone);

fg.Color = bkg;
y=round(sin(az*pi/180)*cos(el*pi/180),2);
x=round(cos(az*pi/180)*cos(el*pi/180),2);
z=round(sin(el*pi/180),2);
xyz = [x,y,z];
dim = find(xyz);
x = xlim;y=ylim;
x = [x(1) (x(2)-x(1))/2 x(2)];y = [y(1) (y(2)-y(1))/2 y(2)];
addnotation=0;
if addnotation
switch dim
    case 1
        if xyz(dim)>0;text([x x(2)],[y(2) y],{'P','S','A','I'},'color',fc,'fontsize',30,'Fontweight','bold','HorizontalAlignment','center')
        else; text([x x(2)],[y(2) y],{'A','S','P','I'},'color',fc,'fontsize',30,'Fontweight','bold','HorizontalAlignment','center');end
    case 2;text([x x(2)],[y(2) y],{'L','S','R','I'},'color',fc,'fontsize',30,'Fontweight','bold','HorizontalAlignment','center')
    case 3;text([x x(2)],[y(2) y],{'L','A','R','P'},'color',fc,'fontsize',30,'Fontweight','bold','HorizontalAlignment','center')
end
end
tmp = getframe(gcf);
img = tmp.cdata;
if nargout>1
imagesc(mask,'parent',ax2);axis(ax2,'off');
fg.Color = [0,0,0];
tmp = getframe(gcf);
mask = rgb2gray(tmp.cdata);
if isfield(im,'mmap')
    cmap= im.mmap;
end
end
end

