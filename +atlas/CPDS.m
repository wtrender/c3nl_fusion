function [L,stats] = CPDS(Y,r,seed,tdarts,ndarts,output)
%% ATLAS.CPDS: Constrained Poisson Disk Sampling for random uniform 
%   segmentation of binary masks
%
%   __           _             
%  / _|         (_)            
% | |_ _   _ ___ _  ___  _ __    
% |  _| | | / __| |/ _ \| `_ \    :- Functional and Structural 
% | | | |_| \__ \ | (_) | | | |      Integration of Neuroimages
% |_|  \__,_|___/_|\___/|_| |_|
%
%
% AUTHOR:  Eyal Soreq
%  EMAIL:  e.soreq14@imperial.ac.uk
%  AFFILIATION:  Imperial College London
%  VERSION:  1.0 CREATED: 03-Nov-2017 08:18:23
%
% INPUTS:
%    Y - some volume to segment (usually 3D mask)
%    r - radius of minimum boundary
%    seed - some number to allow reproducibility
%    dt - dart tries - the maximum number to try to fill in a cell
%    tdarts - dart attempts per cell - the number of attempts before we
%    give up
%    ndarts - number of darts to throw in each itteration
%    output - 'int' (default) or 'binary' - will control the type of output
%              either int which means a integer label per parcel or binary 
%              where the labels are projected into a an additional dimension 
% OUTPUT:
%    L = label map based on output
%    stats = centroid and volume of parcel
%
% DEPENDENCIES:
%   
% This file is part of Fusion Pipeline
% Fusion Pipeline is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% Fusion Pipeline is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with Fusion Pipeline.If not, see <http://www.gnu.org/licenses/>.
%------------- BEGIN CODE --------------
%

if ~exist('seed','var');rng(22041975);else;rng(seed);end
if ~exist('tdarts','var');tdarts = 5;end
if ~exist('ndarts','var');ndarts = 50;end
if ~exist('output','var');output = 'int';end
dim = size(Y);
% pad image to account for bounderies
M = padarray(Y,ones(numel(dim),1)*r,nan); 
Pdim =size(M);
nd = numel(Pdim);
% guess optimal packing for a square to create some upper limit
k = floor(nnz(~isnan(Y))/4*sqrt(2)*r^3);

cellsize = arrayfun(@(x) r+1:round((2*r+1)/sqrt(nd)):x-r-1,Pdim,'un',0);
% generate grid
[grid{1:nd}] = ndgrid(cellsize{:});   
% convert to vector form
sGrid = cellfun(@(x) x(:),grid,'un',0);
sGrid = [cell2mat(sGrid),sub2ind(Pdim,sGrid{:})];

nix = neighb(Pdim,r,'city');
id = size(sGrid);
for ii=1:id(1)
    sGrid(ii,id(2)+1) = nansum(M([sGrid(ii,id(2));sGrid(ii,id(2))+nix]))/numel(nix);
end

ix =find(sGrid(:,id(2)+1)~=0);nnz(ix);
eGrid = false(size(sGrid,1),1); % full Grids
eGrid(ix) = 1; % make only eligible accsessible
ne = nnz(eGrid); % number of Eligible
hGrid = zeros(size(sGrid,1),1); % store dart attempts per cell
np = 0; % number of points
 
pts = [];
ptsix = [];
iter = 0;

tic
rix1 = c3nl.neighb(Pdim,r,'euc');
while np<k && ne > 0
    freeCells = find(eGrid);   %Eligible grids for dart throw  
    d2t = min([ne,ndarts]); % Darts to be thrown
    p = datasample(freeCells,d2t,'Replace',false); % randomly select cells to throw darts in
    tmpPts = cell(numel(p),4);
    for ii=1:numel(p)
        tix = sGrid(p(ii),4)+[0;rix1];
        tix = tix(~isnan(M(tix)));
        tix = tix(randi(numel(tix),1));
        [tmpPts{ii,1:3}] = ind2sub(Pdim,tix);% dart throw on one eligble point from each cell
        tmpPts{ii,4} = tix;
    end
    temp = cell2mat(tmpPts(:,1:3));
    [~,D] = knnsearch([pts;temp],temp,'k',2); %Finding distance between all darts(pts)
    badthrows = D(:,2)<r;
    hGrid(p(badthrows)) = hGrid(p(badthrows))+1; % update misses
    eGrid(p(~badthrows)) = 0; % mark cells as filled
    eGrid = eGrid & (hGrid<tdarts);
    ne = nnz(eGrid);
    pts = [pts;temp(~badthrows,:)]; % add good points to stack
    ptsix = [ptsix;cell2mat(tmpPts(~badthrows,4))];
    np = size(pts,1);
    iter = iter+1;    
    if ~mod(iter,50);fprintf('\n');else;fprintf('*');end
end

L = zeros(Pdim);
C = M;
C(ptsix)=nan;
L(ptsix)=1:numel(ptsix);
for ii=1:r
    lix = neighb(Pdim,ii,'euc');
    for jj=1:numel(ptsix)
        ix = ptsix(jj)+lix;
        ix = ix(~isnan(M(ix(:))));
        ix = ix(L(ix)==0);
        L(ix)= L(ptsix(jj));
        C(ix)=nan;
    end
end


lix = neighb(Pdim,1,'euc');
idx = find(~isnan(C));
for ii=1:numel(idx)
    ix = idx(ii)+lix;
    ix = ix(~isnan(M(ix(:))));
    L(idx(ii)) = mode(L(ix));
end



for ii=1:numel(dim)
    switch ii
        case 1;L = L(1+r:end-r,:,:);
        case 2;L = L(:,1+r:end-r,:);
        case 3;L = L(:,:,1+r:end-r);
    end
end

L = reLabel(L);
stats = regionprops3(L,{'Volume','Centroid'});

switch output
    case 'int';fprintf('\n Output label matrix');
    case 'binary';fprintf('\n Output binary matrix');L = label2Dummy(L);
end

fprintf('\n %d parcels created, with %i mean volume. Total time taken: %0.2f seconds \n ',height(stats),round(mean(stats.Volume)),toc);

end

function [D,W,w]  = neighb(dim,r,dist)

% get back linear shift for neighbourhood

if ~exist('dist','var');dist='euc';end
D = [];
[D{1:length(dim)}]=ndgrid(-r:r);
d = cellfun(@(x) x+1+r,D,'un',0);
d=sub2ind(dim,d{:});
d=d-d((numel(d)+1)/2);
switch dist
    case 'euc';w = (bwdist(d==0,'euclidean'));
    case 'city';w = (bwdist(d==0,'cityblock'));
    case 'chess';w = (bwdist(d==0,'chessboard'));
    case 'qeuc';w = (bwdist(d==0,'quasi-euclidean'));
end

D = d(w<=r*sqrt(2));
W = double(1./w(w<=r*sqrt(2)));W(D==0)=[];
D(D==0)=[];
w = (w>=0)./numel(w);
end


function l = reLabel(L)

idL = unique(L(~isnan(L))); % you atlas id's
idL(idL==0)=[];% ignore zero's
nL = numel(idL);
if max(idL)>=nL
    LUT = zeros(2^16,1,'uint16');
    LUT(idL+1) = 1:nL;
    l = intlut(uint16(L),LUT);
else 
    l=L;
end
end


function [D,nc] = label2Dummy(L)
dim= size(L);
L = double(L(:));
L((isnan(L)))=0;
nc = unique(L(L~=0));
nc(nc==0) = [];
D = false([size(L(:),1),numel(nc)]);
for ii=1:numel(nc)
    D(:,ii) = L==nc(ii);
end
D = reshape(D,[dim,numel(nc)]);

end
%------------- END OF CODE --------------
