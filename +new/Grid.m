classdef Grid
    %% NEW.GRID: creates a new Grid class that stores usefull information about
    %            gridded volumes and adds usefull methods 
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
    %  VERSION:  0.0 CREATED: 23-May-2017 11:42:55
    %
    %% TODO:
    %    1. Extend Grid to hold additinal important info (phase for dwi)
    %    and intensity range for bold
    %    also deal with other volume loaders such as giffi,mif etc..
    %    2. Add Error handling now that this is going to the world
    %   
    %
    %% INPUTS:
    %    Grid(nifti file)
    %
    %% METHODS:
    %  converts from xyz space to native grid
    %  C = convert(G,[1,2,3;5,20,3])
    %
    %% OUTPUT:
    %    Grid - the grid class
    %
    %% EXAMPLES:
    %{
% To get a grid for a nifit file just run:
  G = Grid(nifti(fn));
    %}
    %
    %% DEPENDENCIES: none
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
    properties;
        mm;
        bb;
        mat;
        mat0;
        d;
        vx;
        dimInfo;
        dimOrder
        stride;
        datatype;
    end
    methods
        
        function obj = Grid(G,ext)
            
            switch ext
                case 'nii'
                    mat = G.mat;
                    Q = zeros(4);
                    Q(1:3,4) = -ones(3,1);
                    mat0 = inv(inv(mat)+Q);
                    d = size(G.dat());
                    corners = [ones(1,4),ones(1,4)*d(1);repmat([1,1,d(2),d(2)],1,2);repmat([1,d(3)],1,4);ones(1,8)];
                    XYZ = mat(1:3, :) * corners;
                    bb = [XYZ(:,1),XYZ(:,end)];
                    vx = zeros(3,1);
                    for ii=1:3
                        [~,ix] = max(abs(mat(ii,1:3)));
                        vx(ii) = mat(ii,ix);
                    end
                    obj.mm{1} = linspace(bb(1,1),bb(1,2),d(1));
                    obj.mm{2} = linspace(bb(2,1),bb(2,2),d(2));
                    obj.mm{3} = linspace(bb(3,1),bb(3,2),d(3));
                    obj.bb = bb;
                    obj.mat = mat;
                    obj.mat0= mat0;
                    obj.d = d;
                    obj.vx = vx;
                    obj.dimInfo = G.diminfo;                   
                    [~,or]=max(abs(mat(1:3,1:3)));
                    obj.dimOrder = or;
                    obj.stride = sign([mat(or(1),1),mat(or(2),2),mat(or(3),3)]);
                    obj.datatype = G.dat.dtype;

                case 'mgh'
                    
                case 'mif'
                    mat = G.transform;
                    Q = zeros(4);
                    Q(1:3,4) = -ones(3,1);
                    mat0 = inv(inv(mat)+Q);
                    d = G.dim;
                    corners = [ones(1,4),ones(1,4)*d(1);repmat([1,1,d(2),d(2)],1,2);repmat([1,d(3)],1,4);ones(1,8)];
                    XYZ = mat(1:3, :) * corners;
                    bb = [XYZ(:,1),XYZ(:,end)];
                    vx = G.vox';
                    obj.mm{1} = linspace(bb(1,1),bb(1,2),d(1));
                    obj.mm{2} = linspace(bb(2,1),bb(2,2),d(2));
                    obj.mm{3} = linspace(bb(3,1),bb(3,2),d(3));
                    obj.bb = bb;
                    obj.mat = mat;
                    obj.mat0= mat0;
                    obj.d = d;
                    obj.vx = vx;
                    obj.dimInfo = G.dim	;                   
                    obj.datatype = G.datatype;
                    [~,or]=max(abs(mat(1:3,1:3)));
                    obj.dimOrder = or;
                    obj.stride = G.layout;
            end
        end

        function C = convert(obj,xyz)
            if size(xyz,1)>1
                C = zeros(size(xyz,1),3);
                for ii=1:size(xyz,1)
                    C(ii,:) = [obj.mm{1}(xyz(ii,2)),obj.mm{2}(xyz(ii,1)),obj.mm{3}(xyz(ii,3))];
                end
            else
                C = [obj.mm{1}(xyz(2)),obj.mm{2}(xyz(1)),obj.mm{3}(xyz(3))];
            end
        end
        
    end
end