function im = hemisphere(fc,vx,C,fig,str)
%% PLOT.HEMISPHERE: One line description of what the function or script performs
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
%  VERSION:  0.0 CREATED: 30-Dec-2017 21:03:58
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
hemisphere(input01,input02,input03,input04)
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
disp('')
if ~exist('fig','var');fg = figure(1);else;fg = figure(fig);end
clim = quantile(C,[0.01,0.99]);
C(isnan(C))=-inf;
clf;
x = [0,1.4,4.4,7.4,1.4,4.4]./9+0.025;
y = [1.5,1.5,1.5,1.5,0,0]./5+0.1;
w = [1.4,3,3,1.4,3,3]./9-0.025;
h = [2,2,2,2,1,1]./3.5-0.025;
az = [-180,-90,90,0,-90,90];
el = [0,0,0,0,90,-90];
for ii=1:6
    ax.(['p' num2str(ii)])= axes('Position',[x(ii),y(ii),w(ii),h(ii)]);
    ax.(['p' num2str(ii)]).CLim = clim;
    plot.surface(fc,vx,C,ax.(['p' num2str(ii)]));
    axis tight
    view(az(ii),el(ii));
    camlight;
end
ax.p7 =  axes('Position',[3/9,0.4,0.32,0.03]);
imagesc(linspace(clim(1),clim(2),100),ones(100,1),linspace(clim(1),clim(2),100));
ax.p7.YTick = [];
ax.p7.XTick = linspace(clim(1),clim(2),10);
ax.p7.XTickLabel = round(linspace(clim(1),clim(2),10),1);
title(str,'FontSize',15)
im = getframe(fg);


%------------- END OF CODE --------------
