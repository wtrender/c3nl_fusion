function im = imOverlay(bkim,fgim,Alpha)
%% PLOT.IMOVERLAY: One line description of what the function or script performs
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
%  VERSION:  0.0 CREATED: 18-Aug-2017 19:13:13
%
%% INPUTS:
%    bkim - background image a 2d m x n image with gray values 
%    fgim - a 3d matrix of m x n x k (k can be 1 and then it is just two
%    images overlay)
%    Alpha - the overlay transpancy 
%
%
%% OUTPUT:
%
%% EXAMPLES:
%{
im = plot.imOverlay(bkim,fgim,Alpha)
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
clf;
ax1 = axes('Position',[0 0 1 1]);axis (ax1,'off');
ax2 = axes('Position',[0 0 1 1]);axis (ax2,'off');
imagesc(ax1,bkim);
colormap(ax1,'gray');
axis(ax1,'image');
imagesc(ax2,fgim,'AlphaData',(fgim>0)*Alpha);
%cmap = jet(numel(unique(fgim(:))));
%colormap(ax2,cmap(randperm(size(cmap,1)),:));
axis(ax2,'off')
axis(ax2,'image')
linkaxes([ax1,ax2],'xy');
%------------- END OF CODE --------------
