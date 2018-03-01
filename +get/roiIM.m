function out = roiIM(roi,xyz,dim,lmap,G,ids,cr,bkg)
h = figure(99);clf;
set(h,'units','pixels','position',[0 0 800 800],'Visible','on','resize','off');

ax=axes('position' ,[0 0 1 1]);axis off;
 
switch dim
    case 'a';az = 270;el=90;%axial
    case 'rs';az = 360; el=360;%sagittal
    case 'ls';az = 180; el=360;%sagittal
    case 'lsm';az = 360; el=360;%sagittal
    case 'rsm';az = 180; el=360;%sagittal
    case 'c';az = 270; el=360;%coronal
end

if ~exist('r','var');r=3;end
load +atlas/cortex

switch dim
    case 'rsm';V = cortex.lvol.V;F = cortex.lvol.F;AO = cortex.lvol.AO;
    case 'lsm';V = cortex.rvol.V;F = cortex.rvol.F;AO = cortex.rvol.AO;
    otherwise;V = cortex.V;F = cortex.F;AO = cortex.AO;
end

for jj=1:3
    ax=axes('position',[0,0,1,1]);
    switch jj
            case 1
                patch('Vertices',V,'Faces',F,'FaceColor','w','EdgeColor','none','FaceAlpha',1);
                view(az,el);axis 'off';daspect([1 1 1]);
                lim = axis;
                if strcmp(dim,'lsm');lim(3)=-100;end
                if strcmp(dim,'rsm');lim(4)=100;end
                tmp = getframe(h);
                Output.mask = (rgb2gray(tmp.cdata)~=240).*1;
                cla(ax);
            case 2
                cla(ax);
                p = patch('Vertices',V,'Faces',F,'FaceColor','w','EdgeColor','none','FaceAlpha',0);
                plot.collateSurface(roi,G,0.6,1,lmap,ax);
                ax.YDir='reverse';daspect([1 1 1]);
                view(az,el);axis (ax,lim);axis 'off';
                camlight;lighting phong;
                p.Visible = 'off'
                tmp = getframe(h);
                Output.roi = tmp.cdata;
                Output.alpha = rgb2gray(tmp.cdata)~=240;
                cla(ax);
            case 3; 
                patch('Vertices',V,'Faces',F,'FaceVertexCData',1-AO,'FaceColor','interp','EdgeColor','none','FaceAlpha',1);
                view(az,el);axis 'off';daspect([1 1 1]);colormap gray;axis (ax,lim);
                light
                tmp = getframe(h);
                Output.cortex = tmp.cdata;
                cla(ax);
                patch('Vertices',V,'Faces',F,'FaceVertexCData',AO,'FaceColor','interp','EdgeColor','none','FaceAlpha',1);
                view(az,el);axis 'off';daspect([1 1 1]);colormap gray;axis (ax,lim);
                tmp = getframe(h);
                Output.Ao = rgb2gray(tmp.cdata);
                cla(ax);
    end
end
clf;
if cr
ax0 = axes('position' ,[0 0 1 1]);
imagesc(cat(3,Output.mask.*bkg(1),Output.mask.*bkg(2),Output.mask.*bkg(3)),'AlphaData',Output.mask*0.25,'parent',ax0);
axis(ax0, 'off');axis(ax0,'image');
end
ax1 = axes('position' ,[0 0 1 1]);
imagesc(Output.roi,'AlphaData',Output.alpha,'parent',ax1);
axis(ax1, 'off');axis(ax1,'image');
Alpha = 0.5;
if cr
ax2 = axes('position' ,[0 0 1 1]);
imagesc(Output.cortex,'AlphaData',c3nl.scale(double(Output.Ao)).*Alpha,'parent',ax2);
axis(ax2, 'off');axis(ax2,'image');
colormap(ax2,bone);
end
if ~isempty(ids)
    ax3 = axes('position' ,[0 0 1 1]);
    for ii=1:size(xyz,1)
        text(xyz(ii,2),xyz(ii,1),xyz(ii,3),char(9673),'FontUnits','normalized','FontSize',0.065,'HorizontalAlignment','center','color',lmap(ii,:).^2,'VerticalAlignment','middle')
    end
    axis (ax3,lim); view(az,el);axis 'off';daspect([1 1 1])
    ax3.YDir = 'reverse';
    ax4 = axes('position' ,[0 0 1 1]);   
    for ii=1:size(xyz,1)
        t = rgb2lab(lmap(ii,:));
        if t(1)>70;fontc = 'k';else; fontc = 'w';end
        text(xyz(ii,2),xyz(ii,1),xyz(ii,3),ids(ii),'FontUnits','normalized','FontSize',0.025,'FontWeight','bold','HorizontalAlignment','center','color',fontc,'VerticalAlignment','middle')
    end
    axis (ax4,lim); view(az,el);axis 'off';daspect([1 1 1])
    ax4.YDir = 'reverse';
end


tmp = getframe(h);
out.im = tmp.cdata;
out.mask = rgb2gray(tmp.cdata)~=240;
close( h)
end