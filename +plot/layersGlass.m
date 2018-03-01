function [Out,mmap] = layersGlass(Y,G,dim,thr,cmap,w,h,side,bkg,cr)
%% PLOT.LAYERSGLASS: One line description of what the function or script performs
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
%  VERSION:  0.0 CREATED: 23-Oct-2017 15:22:49
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
layersGlass(input01,input02,input03,input04)
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
d = size(Y);
if ~exist('side','var')||isempty(side);side = 3;end
if ~exist('bkg','var')||isempty(bkg);bkg=[236,246,256]./256;end

[im,mmap] = plot.mip(Y,G,dim,thr,cmap); 
if ~exist('fig','var');fig = 1;end
fg = figure(fig);clf;
set(fg,'units','pixels','position',[0 0 w h],'Visible','on','resize','off');
if ~cr
ax0 = axes('position' ,[0.02 0.05 0.9 0.9]);
imagesc(im.mask,'AlphaData',im.mask*0.25,'parent',ax0);
axis(ax0, 'off');axis(ax0,'image');
Alpha = 0.5;
colormap(ax0,bkg);
end
ax1 = axes('position' ,[0.02 0.05 0.9 0.9]);
imagesc(im.img,'AlphaData',im.alpha.*im.mask,'parent',ax1);
axis(ax1, 'off');axis(ax1,'image');

if ~cr
ax2 = axes('position' ,[0.02 0.05 0.9 0.9]);
imagesc(im.cortex,'AlphaData',c3nl.scale(double(im.Ao)).*im.mask.*Alpha,'parent',ax2);
axis(ax2, 'off');axis(ax2,'image')
colormap(ax2,bone);
linkaxes([ax0,ax1,ax2])
end

xy = [find(any(im.mask,1),1,'first'),find(any(im.mask,1),1,'last'),find(any(im.mask,2),1,'first'),find(any(im.mask,2),1,'last')];

axis(xy)
tmp = getframe(fg);
Out.im = tmp.cdata;
Out.mask = rgb2gray(tmp.cdata)~=240;

end



%------------- END OF CODE --------------
