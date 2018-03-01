function ax = mypcshow(ax,vol,G,thr,marker,cmap,bkg)
% mypcshow([],pos_T,G,4,[],[],hot)
if ~exist('marker','var')||isempty(marker);marker='.';end
step = 1;
if isempty(G)
   d = size(vol); 
   G = struct('xmm',1:d(1),'ymm',1:d(2),'zmm',1:d(3)); 
end
if ~exist('bkg','var')||isempty(bkg)
    bkg = [0,0,0];
end
[x,y,z] = meshgrid(G.mm{2},G.mm{1},G.mm{3});

ix = vol <= thr | isnan(vol);

I = vol(~ix);
% II = vol;
% II(ix) = NaN;
if isempty(ax);ax = newplot;end
% h = slice(ax,x,y,z,II,min(G.ymm):step:max(G.ymm),min(G.xmm):step:max(G.xmm),min(G.zmm):step:max(G.zmm));
% for ii=1:size(h,1)
%     h(ii).FaceColor = 'interp';
%     h(ii).EdgeColor = 'none';
% end
scatter3(ax,x(~ix)-randn(sum(~ix(:)),1),y(~ix)-randn(sum(~ix(:)),1),z(~ix)-randn(sum(~ix(:)),1),'sizeData',round(I.^3)./2,'CData',I,'marker',marker)
axis 'off'
hFigure = get(ax,'Parent');
vertAxis = 'Z';
vertAxisDir = 'Up';
% Lower and upper limit of auto downsampling.
ptCloudThreshold = [1920*1080, 1e8]; 

% Initialize point cloud viewer controls.
vision.internal.pc.initializePCSceneControl(hFigure, ax, vertAxis,...
    vertAxisDir, ptCloudThreshold, true);
colormap(cmap);
alpha(0.5);
hFigure.Color = bkg;
view(-90,90)
end