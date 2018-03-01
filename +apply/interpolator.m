function [V,G] = interpolator(source,target,type,output)
%% APPLY.INTERPOLATOR: One line description of what the function or script performs
%
%    __           _             
%   / _|         (_)            
%  | |_ _   _ ___ _  ___  _ __    
%  |  _| | | / __| |/ _ \| '_ \    :- Functional and Structural 
%  | | | |_| \__ \ | (_) | | | |      Integration of Neuroimages
%  |_|  \__,_|___/_|\___/|_| |_|
%
%
%% AUTHOR:  Eyal Soreq
%  EMAIL:  e.soreq14@imperial.ac.uk
%  AFFILIATION:  Imperial College London
%  VERSION:  0.0 CREATED: 24-May-2017 06:02:24
%

%% TODO: 
%   1. add save.vol 
%

%% INPUTS:
%    source - 
%    target - 
%    type - 
%    output - 
%
%
%% OUTPUT:
%    V - 
%    G - 
%
%% EXAMPLES:
%{
[V,G] = interpolator(source,target,type,output)
%}
%
%% DEPENDENCIES: save.vol, load.vol
%
% This file is part of C^3NL Pipeline
% C^3NL Pipeline is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% C^3NL Pipeline is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with C^3NL Pipeline.If not, see <http://www.gnu.org/licenses/>.
%------------- BEGIN CODE --------------
%
if ~exist('type','var');type='nearest';end
st = find([ischar(source),isstruct(source),iscellstr(source)]);
if st==1;[sv,sg] = load.vol(source);end
if st==2;sv=source.V;sg=source.G;end
if st==3;sv=source{1};sg=source{2};end
tt = find([ischar(target),isstruct(target),iscellstr(target)]);
if tt==1;[tv,tg] = load.vol(target);end
if tt==2;tv=target.V;tg=target.G;end
if tt==3;tv=target{1};sg=target{2};end
clear tt st
%% get size of source & target
sd = sg.d;
td = tg.d(1:3);

%% get dimensional order of source and target
smat = sg.mat;
tmat = tg.mat;
if ~exist('spm_matrix')
s = strsplit(which(['apply.' mfilename]),'+');
addpath([s{1} '3rdParty' filesep 'spm12']);
end

%% apply 3D interpolation based on type
switch type
    case 'reslice'
        V = zeros(tg.d(1:3));
        for ii = 1:tg.d(3)
            M = inv(spm_matrix([0 0 -ii])*inv(tmat)*smat);
            V(:,:,ii) = spm_slice_vol(sv, M, tg.d(1:2), 0); %
        end
        G=tg;
    case 'reslice_roi'
        nc = unique(sv(~isnan(sv)));nc(nc==0) = [];
        V = zeros([tg.d(1:3),numel(nc)]);
        for ii=1:numel(nc)
            D1 = sv==nc(ii);
            l1 = apply.interpolator(struct('V',double(D1),'G',sg),struct('V',tv,'G',tg),'reslice');
            ix=find(l1);
            [~,~,K]=ind2sub(size(l1),ix);
            for z=min(K):max(K)                    
                M = inv(spm_matrix([0 0 -z])*inv(tmat)*smat);
                V(:,:,z,ii) = spm_slice_vol(double(D1), M, tg.d(1:2), 1); %
            end
            if mod(ii,80);fprintf('%%');else;fprintf('-%i\n',ii);end
        end
        G=tg;
        G.d = size(V);
    case {'nearest', '*linear','linear','cubic','*cubic','spline'}
        [~,sdim]=max(abs(smat(1:3,1:3)));
        [~,tdim]=max(abs(tmat(1:3,1:3)));
        %% find dimension to flip to ensure monotonic increasing grids in both target and source
        mat =smat(1:3,1:3);
        sflip = find(mat(sub2ind([3,3],1:3,sdim))<0);
        mat =tmat(1:3,1:3);
        tflip = find(mat(sub2ind([3,3],1:3,tdim))<0);

        %% flip dimesions if needed 
        [fsv,fsg]= apply.flipVol(sv,sg,sflip);
        [~,ftg]= apply.flipVol(tv,tg,tflip);

        [Xs,Ys,Zs] = meshgrid(fsg.mm{2},fsg.mm{1},fsg.mm{3});
        [Xt,Yt,Zt] = meshgrid(ftg.mm{2},ftg.mm{1},ftg.mm{3});       
        
        
        if numel(sd)==3
            V = interp3(Xs,Ys,Zs,fsv,Xt,Yt,Zt,type);           
            G = ftg;
            if ~isempty(tflip)
                [V,G]= apply.flipVol(V,ftg,tflip); % flip back to match target grid
            end
        else
            V = zeros([td(1:3),sd(4)]);
            for ii=1:sd(4)
                V(:,:,:,ii) = interp3(Xs,Ys,Zs,fsv(:,:,:,ii),Xt,Yt,Zt,type);           
            end
            G = ftg;
            G.d = size(V);
        end
end
%rmpath([source{1} '3rdParty' filesep 'spm12'])

if exist('output','var')
   fast_nii_save(V,G,output); 
end

end



%------------- END OF CODE --------------
