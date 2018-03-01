function adjacency(fs,cmap,T,or,outdir,fn,w,h,thr)
disp('')
if ~exist('thr','var');thr=0;end
%{
mdl = tmp.data.mdl;
nnz(mdl.beta{3}(:,30))
%mdl =data.mdl;
cmap = [151,84,204;129,198,121;243,184,99]./256;% frac/Num/Pos

tmap = hsv(360);
fmap = genCmap([tmap(220,:);tmap(220,:).^.5;.85,.85,.85;cmap(1,:);cmap(1,:).^2]);
nmap = genCmap([tmap(220,:);tmap(220,:).^.5;.85,.85,.85;cmap(2,:);cmap(2,:).^3]);
pmap = genCmap([tmap(220,:);tmap(220,:).^.5;.85,.85,.85;cmap(3,:);cmap(3,:).^3]);

fs = mdl.beta{3};
plot.adjacency(fs,pmap,T,'a',outdir,'POS');
plot.adjacency(fs,pmap,T,'rs',outdir,'POS');
plot.adjacency(fs,pmap,T,'ls',outdir,'POS');
fs = mdl.beta{2};
plot.adjacency(fs,nmap,T,'a',outdir,'NUM');
plot.adjacency(fs,nmap,T,'rs',outdir,'NUM');
plot.adjacency(fs,nmap,T,'ls',outdir,'NUM');

fs = mdl.beta{1};
plot.adjacency(fs,fmap,T,'a',outdir,'FRAC');
plot.adjacency(fs,fmap,T,'rs',outdir,'FRAC');
plot.adjacency(fs,fmap,T,'ls',outdir,'FRAC');
INPUT
=====
take adj matrix, 3d grid, colormap, information table per nodes, 3d
orientation and plot the 3d graph. 

OUTPUT
======
im = 1 or more (animation) image in the specific orientation of the
neuronal graph

%}
%adjPlot
bkg=[230,250,250]./256;
load atlas/layers96
ld = layers.(or);
%mask = ld.Alpha;
az = ld.az;
el = ld.el;
side = T.regionCode;
side(side==0)=nan;
side = mod(side,2);
mask = side + side';
m = height(T);
f =  figure(999);clf
f.Color = [0,1,0];
for ii=1:3
    id = ['p' num2str(ii)];
    ax.(id) = axes('position',[0.05 0.05 0.9 0.9]);axis(ax.(id), 'off');
end

imagesc(ld.Alpha,'AlphaData',ld.Alpha,'parent',ax.p1);
colormap(ax.p1,bkg);axis(ax.p1, 'image');axis(ax.p1, 'off');
imagesc(ld.ao.img,'AlphaData',(ld.ao_alpha.img.*ld.Alpha*.65),'parent',ax.p3);
colormap(ax.p3,bone);axis(ax.p3, 'off');
axis(ax.p3, 'image');
for ii=1:size(fs,2)
    switch or
        case 'a'; Output = genImage(fs(:,ii),T,az,el,true(m),cmap,thr);                   
        %case 'rs';Output = genImage(fs(:,ii),T,az,el,mask==0,cmap);  
        case 'rs';Output = genImage(fs(:,ii),T,az,el,true(m),cmap,thr);  
        %case 'ls';Output = genImage(fs(:,ii),T,az,el,mask==2,cmap); 
        case 'ls';Output = genImage(fs(:,ii),T,az,el,true(m),cmap,thr);  

    end
    cla(ax.p2);
    if ~isempty(Output)
        imagesc(Output.img,'AlphaData',Output.alpha,'parent',ax.p2);
        axis(ax.p2, 'image');axis(ax.p2, 'off');
    end
    figure(999);
    save.pdf(sprintf('%s%s%s_%02d_%s',outdir,filesep,or,ii,fn),w,h)
    %print('-dpng',sprintf('%s%s%s_%s_%02d.png',outdir,filesep,fn,or, ii));
end

%ffmpegimages2video([outdir,filesep,fn,'_',or,'_','%02d.png'],[outdir,filesep,fn,'_',or,'.gif'],'InputFrameRate',5,...
%      'VideoCodec','gif','DeleteSource','off');
end


function Output = genImage(fs,T,az,el,mask,cmap,thr)
m = height(T);
ix=find(triu(true(m),1));
adj = zeros(m);
adj(ix(fs~=0))=fs(fs~=0);
adj = symAdj(adj,'upper');
nw = sum(adj)';
if nnz(nw~=0)
clim = round(max(abs(nw))).*[-1,1];
nc = ind2Cmap(cmap,[clim(:);nw]);
T.nmap = nc(3:end,:);
adj(adj<thr)= 0;
adj = adj.*mask;
Output = roi2imag([],struct('type','adj','adj',adj,'az',az,'el',el,'dpi',96,'T',T,'cmap',cmap));
else 
    Output =[];
end
end