function splines2d(X,Y,ax,cmap,labels,E,ms,lim,fs,tr)
%% PLOT.SPLINES2D: plot a curved plot with or without errors 
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
%  VERSION:  0.0 CREATED: 25-Oct-2017 07:50:56
%
%% INPUTS:
%    X - spacing  
%    Y - magnitude
%    ax - axis object 
%    cmap - color per row 
%    labels - axis labels    
%    E - error bounderies   
%% OUTPUT:
%
%% EXAMPLES:
%{
plot.splines2d(X,Y,ax,cmap,labels,E)
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
if ~exist('ax','var');ax=gca;end
if ~exist('cmap','var');cmap=jet;end

if ~exist('ms','var');ms=10;end
if ~exist('tr','var');tr=0.2;end
hold(ax, 'on');
d= size(Y);
xx = X(1):min(diff(X))/10:X(end);
for ii=1:d(1)
    yy = pchip(X,Y(ii,:),xx);
    %yy = spline(X,Y(ii,:),xx);
    if exist('E','var')
        yyp = pchip(X,Y(ii,:)+E(ii,:),xx);
        yym = pchip(X,Y(ii,:)-E(ii,:),xx);
        patch([xx,fliplr(xx)],[yyp,fliplr(yym)],cmap(ii,:),'FaceAlpha',tr,'EdgeColor','none','parent',ax)
    end
    plot(xx,yy,'Color',cmap(ii,:).^4,'parent',ax);   
    plot(X,Y(ii,:),'o','MarkerFaceColor',cmap(ii,:),'MarkerSize',ms,'MarkerEdgeColor','none','parent',ax);
     
end
axis tight;
if ~exist('lim','var')||isempty(lim)
    if exist('E','var')
        lim = [X(1)-min(diff(X))/10,X(end)+min(diff(X))/10,...
               min(Y(:))-abs(max(E(:))),max(Y(:))+abs(max(E(:)))];
        
    else
        lim = [X(1)-min(diff(X))/10,X(end)+min(diff(X))/10,...
               min(Y(:))-abs(min(Y(:))*0.1),max(Y(:))+abs(min(Y(:))*0.1)];
    end
end
lim = round(lim,2);
axis(lim);
yy = ylim;
ax.YTick = [yy(1),yy(2)];
ax.XTick = X;
ax.XTickLabel = labels;
ax.FontUnits = 'normalized';
ax.FontSize = fs;
end
%------------- END OF CODE --------------
