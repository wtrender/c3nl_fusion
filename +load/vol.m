function [varargout] = vol(fn,hdronly,ext)
%% LOAD.VOL: detects the type of volume and loads it with the right loader
%    __           _             
%   / _|         (_)            
%  | |_ _   _ ___ _  ___  _ __    
%  |  _| | | / __| |/ _ \| '_ \    :- Functional and Structural 
%  | | | |_| \__ \ | (_) | | | |      Integration of Neuroimages
%  |_|  \__,_|___/_|\___/|_| |_|
%
%% AUTHOR:  Eyal Soreq
%  EMAIL:  e.soreq14@imperial.ac.uk
%  AFFILIATION:  Imperial College London
%  VERSION:  0.0 CREATED: 23-May-2017 18:10:47
%
%% TODO:
%  add error catching 

%% INPUTS:
%    fn - either a char vector with one file name pointing to a 3D/4D volume 
%         or a cell of strings each pointing to a volume 
%         or an SPM12 struct 
%    hdronly - if true only G will be extracted
%           
%% OUTPUT:
%    V - either a 3D/4D matrix or a cell array of matrices 
%    G - either one grid class or a cell array of grids 
%
%% EXAMPLES:
%{
[V,G] = load.vol(fn)
%}
%
%% DEPENDENCIES:
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
if ~exist('hdronly','var');hdronly=false;end
vartype = find([ischar(fn),isstruct(fn),iscellstr(fn)]);
switch vartype
    case 1 % simple file deal with ext
        if ~exist('ext','var');[~,~,ext]=fileparts(fn(1,:));end
        switch ext
            case '.gz'
                system(['gunzip -c ' fn ' > ' fn(1:end-3)]);
                if hdronly;varargout{1} = load.vol(fn(1:end-3),hdronly);
                else;[varargout{1},varargout{2}] = load.vol(fn(1:end-3),hdronly);
                end                
                system(['rm -rf ' fn(1:end-3)]);
            case '.nii'
                source = strsplit(which(['load.' mfilename]),'+');
                addpath([source{1} '3rdParty' filesep 'spm12']);
                NI = nifti(fn);
                if hdronly;varargout{1} =new.Grid(NI(1),'nii');
                else                     
                    varargout{1} = NI.dat();
                    varargout{2} = new.Grid(NI(1),'nii');
                end
                rmpath([source{1} '3rdParty' filesep 'spm12'])
            case '.vtk'
                varargout{1} = load.vtk(fn);
            case 'cifti'
                source = strsplit(which(['load.' mfilename]),'+');
                addpath([source{1} '3rdParty' filesep 'cifti-matlab']);
                varargout{1} = ft_read_cifti(fn);
                
            case '.mif'
                source = strsplit(which(['load.' mfilename]),'+');
                addpath([source{1} '3rdParty' filesep 'mrtrix']);
                tmp = read_mrtrix(fn);
                varargout{1} = tmp.data;
                tmp.data = [];
                varargout{2} = new.Grid(tmp,'mif'); 
                rmpath([source{1} '3rdParty' filesep 'mrtrix'])

            case '.mgz'
                disp('not implemented yet')
                return
            case '.gii'
                source = strsplit(which(['load.' mfilename]),'+');
                addpath([source{1} '3rdParty' filesep 'spm12']);
                NI = gifti(fn);
                varargout{1} = NI;

            
            case '.dcm'
                disp('not implemented yet')
                return    
            case '.IMA'
                disp('not implemented yet')
                return    
    
        end
    case 2
        % assume SPM structure
        % assume identical grid 
        P = fn.xY.P;
        NI = nifti(P);
        V = [];
        if ~hdronly;
            for ii=1:numel(NI)
                tmp = NI(ii).dat();
                V(:,:,:,ii) = tmp;
            end
        end
        G = new.Grid(NI(1),'nii');
    case 3       
        V = cell(numel(fn),1);
        for ii=1:numel(fn) 
            [varargout{1}{ii},varargout{2}{ii}] = load.vol(fn{ii},hdronly);
        end                
    otherwise
        fprintf('%s\n%s\n%s\n','INPUT must be either a char vector with one file name', 'pointing to a 3D/4D volume', 'or a cell of strings each pointing to a volume','or an SPM12 struct')
end

end

%------------- END OF CODE --------------
