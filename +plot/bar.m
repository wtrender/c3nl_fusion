function bar(X,gp,labels,cmap,fig,lim,se)
%% PLOT.BAR: One line description of what the function or script performs
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
%  VERSION:  0.0 CREATED: 20-Jul-2017 12:07:07
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
bar(input01,input02,input03,input04)
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

d = size(X);
if numel(d)<3;d(3) =1;end
if ~exist('gap','var');gap = 0.2;end
if ~exist('cmap','var');cmap = jet(d(2)*d(3));end
if ~exist('fuzzy','var');fuzzy = 0;end
if ~exist('lim','var')||isempty(lim);lim = [0,max(X(:))*1.005];end
if ~exist('yx','var');yx = 0;end
if ~exist('fig','var');figure(1);else ;figure(fig);end
if ~exist('fc','var');fc=0.02;end
if ~exist('se','var');se=false;end

[y,x,h,w] = get.fancyGrid(ones(d(2),1),0.01,0.2,'stretch');

clf;
ax = [];
for ii=1:d(2)
    ax.(['p' num2str(ii)]) = axes('position',[x(ii),y(ii),w,h]);
    grid on    
end
GP = unique(gp);
hold all;
c = 0;
for ii=1:d(2)   
    for jj=1:numel(GP)
        c = c+1;
        ix = gp ==GP(jj);
        xt = X(ix,ii);
        xm = nanmean(xt);
        sx = nanstd(xt);
        if se
           sx = sx./sqrt(se);
        end
        [xx,yy]=plot.square(jj,xm/2,0.5,xm);
        patch(xx,yy,cmap(jj,:),'FaceAlpha',0.75,'parent',ax.(['p' num2str(ii)]),'Edgecolor','k','linewidth',1);
        line([jj,jj],[xm-sx,xm+sx],'linewidth',1,'color','k','parent',ax.(['p' num2str(ii)]));
        line([jj-0.125,jj+0.125],[xm-sx,xm-sx],'linewidth',1,'color','k','parent',ax.(['p' num2str(ii)]));
        line([jj-0.125,jj+0.125],[xm+sx,xm+sx],'linewidth',1,'color','k','parent',ax.(['p' num2str(ii)]));
    end
    line([0,numel(GP)+1],[0,0],'linewidth',0.5,'color','k','parent',ax.(['p' num2str(ii)]),'LineStyle','--');
end

for ii=1:d(2)
    axis(ax.(['p' num2str(ii)]),[0.25,numel(GP)+0.75,lim])
    ax.(['p' num2str(ii)]).YTick =round(linspace(lim(1),lim(2),5),2);
    if ii>1
        ax.(['p' num2str(ii)]).YTickLabel = '';
    end
    ax.(['p' num2str(ii)]).XTick = 1:numel(GP);
    ax.(['p' num2str(ii)]).XTickLabel = char(GP);
    ax.(['p' num2str(ii)]).XTickLabelRotation = -30;
    ax.(['p' num2str(ii)]).Color = [0.9,0.9,0.9]; 
    title(ax.(['p' num2str(ii)]),labels{ii})
end

%------------- END OF CODE --------------
