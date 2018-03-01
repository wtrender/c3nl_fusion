function h = vol(V,varargin)
%% PLOT.VOL: One line description of what the function or script performs
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
%  VERSION:  0.0 CREATED: 20-Aug-2017 15:40:12
%
% INPUTS:
%    V - volume to plot 
%    varargin - paired argument of the following inputs  
%       type - 
%           'mip' = maximum intensity projection
%           'iso' = convert a volume of labels to isosurfeces
%           'pc'  = clouds of points   
%       cmap  - 
%           a matrix in the form of m x 3 RGB values m being the number of
%           labels in the matrix in the case of 'iso' (default = jet)
%           for both mip or pc a gradient is supplied 
%       Clim - color limit in the form of [minimum, maximum] threshold 
%              default = 
%                   +/- maps = max(abs(V(:))*[-1,1] 
%                   otherwise  = max(V(:))*[0,1]
%       grid = a grid class object that matches the volume in question
%       Alpha = A scalar value between 0-1 that controles isosurface
%       transparency 
%       smooth = a scalar controling a 3d gaussian smmothing factor    
%       
%               
%
%
% OUTPUT:
%       h = figure handle     
%
% EXAMPLES:
%{
V =  round(imgaussfilt3(rand(50,50,50)*100));
h = plot.vol(V,'cmap',jet(100))
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
collateSurface
%------------- END OF CODE --------------
