function roi(L,G,T,W,cmap,nodes,cr)
%% PLOT.ROI: One line description of what the function or script performs
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
%  VERSION:  0.0 CREATED: 18-Jun-2017 11:06:07
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
plot.roi(input01,input02,input03,input04)
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
if ~exist('cr','var');cr =1;end
bkg = [0.9,0.9,0.95];
n = height(T);
[y,x,h,w] = get.fancyGrid([2;1;2],0.015,0.025,'weighted',[.3,.45,.3]);
d = size(L);
rs = L;rs(round(d(1)/2)-2:end,:,:,:) = 0;
ls = L;ls(1:round(d(1)/2)+2,:,:,:) = 0;
or = {'ls','rsm','a','rs','lsm'};
load +atlas/cortex
num = arrayfun(@(a) sprintf('%02d',a),T.ROIid,'un',0);
str = cellfun(@(a) sprintf('    %s',a), T.name_AAL2,'un',0);
cmap2 = zeros(n,3);
for ii=1:n
    cmap2(ii,:) = cmap(W(ii),:);
end

im = [];m=[];
for ii=1:numel(or)
    switch or{ii}
        case 'a';roi = L;
        case 'rs';roi = rs;
        case 'ls';roi = ls;
        case 'rsm';roi = rs;
        case 'lsm';roi = ls;
    end
    if numel(d)==3
        ix1=ismember(T.ROIid,double(unique(roi(:))));
    else
        l = reshape(roi,[],d(4));
        ix1 = sum(l)>0;
    end
    N = num(ix1);
    lmap=cmap2(ix1,:);
    xyz = [T.X(ix1),T.Y(ix1),T.Z(ix1)];
    if ~nodes;im.(or{ii}) = get.roiIM(roi,xyz,or{ii},lmap,G,[],cr,bkg);
    else; im.(or{ii}) = get.roiIM(roi,xyz,or{ii},lmap,G,num2cell(T.ROIid(ix1)),cr,bkg);
    end 
end

fg = figure(100);clf
for ii=1:numel(or)
    axh.(['P' num2str(ii)]) = axes('position',[x(ii),y(ii),w(ii),h(ii)]);
    imagesc(im.(or{ii}).im,'AlphaData',im.(or{ii}).mask,'parent',axh.(['P' num2str(ii)]));axis(axh.(['P' num2str(ii)]),'image');
    axis(axh.(['P' num2str(ii)]),[find(any(im.(or{ii}).mask,1),1,'first'),find(any(im.(or{ii}).mask,1),1,'last'),find(any(im.(or{ii}).mask,2),1,'first'),find(any(im.(or{ii}).mask,2),1,'last')]);
    axis off
end
    
    
 

end






%------------- END OF CODE --------------
