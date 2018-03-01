
cd /home/eyalsoreq/Dropbox/PhD/CURRENT/papers/Watershed/Examples/surface
sp = load.vol('hcp.embed.dscalar.nii',0, 'cifti');
Lg = load.vol('S900.L.midthickness_MSMAll.32k_fs_LR.surf.gii');
Rg = load.vol('S900.R.midthickness_MSMAll.32k_fs_LR.surf.gii');

sv = sp.x100307_tfMRI_MOTOR_level2_CUE_hp200_s2;
sv2 = sp.x100307_tfMRI_MOTOR_level2_AVG_hp200_s2;
Lv= sv(1:size(Lg.vertices,1));
Rv= sv(size(Lg.vertices,1)+1:size(Lg.vertices,1)*2);

Lv2= sv2(1:size(Lg.vertices,1));
Rv2= sv2(size(Lg.vertices,1)+1:size(Lg.vertices,1)*2);
%% run in free surfer 
LA = dlmread('L.adj.txt');
RA = dlmread('R.adj.txt');
Al = false(size(LA,1));
Ar = false(size(LA,1));

for ii=1:size(LA,1)
    if LA(ii,2)==5;edges = LA(ii,3:7)+1;else;edges = LA(ii,3:8)+1;end
    for jj=1:numel(edges);Al(LA(ii,1)+1,edges) = 1;end
    if RA(ii,2)==5;edges = RA(ii,3:7)+1;else;edges = RA(ii,3:8)+1;end   
    for jj=1:numel(edges);Ar(RA(ii,1)+1,edges) = 1;end
end
G= graph(Al);
tic;nodeIDs = nearest(G,10,3);toc
tic;id = unique(atlas.adjDLS(Al,10,2));toc
unique(id)
figure(2);imagesc(Al)

Lvc = struct('cdata', Lv);
tic;
Lmx = atlas.AdjWS(Lv2,Al,30,1,'max');
Rmx = atlas.AdjWS(Rv,Ar,20,1,'max');
toc
Lvc = struct('cdata', Lmx);
figure; plot(Lg,Lvc);
colormap(jet)
Rvc = struct('cdata', Rmx);
plot(Rg,Rvc);

cmap = jet(numel(unique(Lmx)));

cmap = cmap(randperm(numel(unique(Lmx))),:)
colormap(cmap)