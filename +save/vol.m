function vol(V,G,fn,ftype,datatype)
%% SAVE.VOL: One line description of what the function or script performs
%
%   __           _             
%  / _|         (_)            
% | |_ _   _ ___ _  ___  _ __    
% |  _| | | / __| |/ _ \| `_ \    :- Functional and Structural 
% | | | |_| \__ \ | (_) | | | |      Integration of Neuroimages
% |_|  \__,_|___/_|\___/|_| |_|
%
%
%% AUTHOR:  Eyal Soreq
%  EMAIL:  e.soreq14@imperial.ac.uk
%  AFFILIATION:  Imperial College London
%  VERSION:  0.0 CREATED: 16-Jun-2017 13:50:14
%
%% INPUTS:
%    input01 - 
%    input02 - 
%    input03 - 
%    input04 - 
%
%
%% OUTPUT:
%
%% EXAMPLES:
%{
vol(input01,input02,input03,input04)
%}
%
%% DEPENDENCIES:
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
% if ~isfield(G,'datatype')
%     datatype = 'float32-le';
% else 
%     datatype = G.datatype;
% end
switch ftype
    case 'mif'
        source = strsplit(which(['save.' mfilename]),'+');
        addpath([source{1} '3rdParty' filesep 'mrtrix']);
        im.data = V; % a N-dimensional array (N <= 16)
        im.vox = G.vx; % N-vector of voxel sizes (in mm) (default: { 2 }) [optional]
        im.comments = {'5tt Exported by Eyal Soreq'};
        im.datatype = datatype; % the datatype specifier (default: float32) [optional]
        im.transform = G.mat; % a 4x4 matrix
        write_mrtrix(im,fn);       
    case 'nii'
        source = strsplit(which(['save.' mfilename]),'+');
        addpath([source{1} '3rdParty' filesep 'spm12']);
        d = size(V);
        dat = file_array(fn,d,datatype,ceil(348/8)*8);
        NO=nifti;
        NO.mat = G.mat;
        NO.mat0 = G.mat0;
        NO.mat_intent = 'Aligned';
        NO.extras = [];
        NO.diminfo = '';
        NO.dat = dat;
        NO.descrip = '';
        create(NO);
        switch numel(d)
            case 2;	NO.dat(:,:) = V;
            case 3;	NO.dat(:,:,:) = V;
            case 4;	NO.dat(:,:,:,:) = V;
            case 5;	NO.dat(:,:,:,:,:) = V;
        end
        rmpath([source{1} '3rdParty' filesep 'spm12'])
    case 'mat'    
    
    
    
    
end

end
%------------- END OF CODE --------------
