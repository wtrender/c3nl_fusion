function fv = collateSurface(L,G,A,s,cmap,ax,clip)

% Input
% =====
% L=binary volume
% G = grid object
% alpha = level of transparency
% smoothness of volume 
% color map 
% figure handle 
% clipping attributes [x,y,z]
rnd =0;
cc = unique(L(L>0));
if ~exist('ax','var')||isempty(ax);ax = gca;end
if ~exist('clip','var')||isempty(clip);clip = 0;end
if nargin<7;clip = 0;end
if size(L,4)>1
    sz = size(L);
    cc = 1:sz(4);
    d4 = 1;
else 
    d4 = 0;
end
if ~exist('A','var');A=ones(numel(cc),1)*0.5;end
if numel(A)~=numel(cc);A = A*ones(numel(cc),1);end
if ~exist('s','var');s=1;end
if ~exist('cmap','var');cmap=jet(numel(cc));end
if ~exist('G','var')||isempty(G)
   d = size(L); 
   G.mm{1} = 1:d(1);
   G.mm{2} = 1:d(2); 
   G.mm{3} = 1:d(3); 
end
if isstruct(clip);ymm = clip.ymm;xmm = clip.xmm;zmm = clip.zmm;end
if ~exist('xmm','var');xmm = G.mm{1};end
if ~exist('ymm','var');ymm = G.mm{2};end
if ~exist('zmm','var');zmm = G.mm{3};end
[X,Y,Z]= meshgrid(ymm,xmm,zmm);
if isempty(cmap);cmap = hot(numel(cc));end
if rnd;cmap = cmap(randperm(numel(cc)),:);end
%smooth = 1;
%debug =1;
axis([min(ymm) max(ymm) min(xmm) max(xmm) min(zmm) max(zmm)])
set(gca,'Ydir','reverse')
for ii=1:numel(cc)
    if ~d4;tmp = L == cc(ii);
    else tmp = L(:,:,:,ii);end 
    if s>0
        tmp = smooth3(tmp,'gaussian',3,s);
    end
    fv = isosurface(X,Y,Z,tmp,.5);
    patch(fv,'FaceColor',cmap(ii,:),'EdgeColor','none',...
        'facealpha',A(ii),'facelighting','gouraud',...
        'SpecularColorReflectance',0.2 ,'DiffuseStrength',0.8,...
        'SpecularExponent',3,'SpecularStrength',0.4,'Parent', ax);
    if exist('debug','var')
        view(-90,0)
        daspect([1 1 1])
        drawnow();
        pause();
    end
end
view(3)
daspect([1 1 1])
axis([min(ymm) max(ymm) min(xmm) max(xmm) min(zmm) max(zmm)]);
ax.YDir = 'reverse';

end