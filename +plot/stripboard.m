function stripboard(X,labels,fig,lim,pl)
%% PLOT.STRIPBOARD: One line description of what the function or script performs
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
%  VERSION:  0.0 CREATED: 21-Jul-2017 04:35:58
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
stripboard(input01,input02,input03,input04)
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

d = size(X);
if numel(d)<3;d(3) =1;end
if ~exist('gap','var');gap = 0.2;end
if ~exist('cmap','var');cmap = jet(d(2)*d(3));end
if ~exist('fuzzy','var');fuzzy = 0;end
if ~exist('lim','var');lim = [0,max(X(:))*1.005];end
if ~exist('yx','var');yx = 0;end
if ~exist('fig','var');figure(1);else ;figure(fig);end
if ~exist('fc','var');fc=0.02;end
if ~exist('pl','var');pl=[1,1,1];end

gmap = flipud(gray);


sb = cat(3,squeeze(mean(X)),squeeze(std(X)));
tmap =  c3nl.ind2Cmap(gmap,[lim(:);sb(:)]);
tmap = reshape(tmap(3:end,:),[size(sb),3]);
ssb = sb;
ssb(:,:,2) = sb(:,:,1)+sb(:,:,2);
ssb = c3nl.scale([lim(:);ssb(:)]);
ssb = reshape(ssb(3:end),size(sb));

clf;
ax =  axes('position',[0.2,0.2,0.6,0.6]);
hold all;
for ii=1:d(2)   
    for jj=1:d(3)
        [x1,y1]=plot.circle(jj,ii,ssb(ii,jj,1)/2);
        [x2,y2]=plot.circle(jj,ii,ssb(ii,jj,2)/2);
        patch(x2,y2,squeeze(tmap(ii,jj,2,:))','edgecolor','k');
        patch(x1,y1,squeeze(tmap(ii,jj,1,:))','edgecolor','w');
        text(jj,ii,num2str(round(sb(ii,jj,1),1)),'FontUnits','normalized','fontsize',0.2/d(2),...
            'HorizontalAlignment','center','VerticalAlignment','middle','color','w');
    end
end
axis(ax,'image');
axis(ax,[0.5,d(3)+0.5,0.5,d(2)+0.5]);
ax.Color = [0.9,0.9,0.9]; 
ax.XTick = 1:d(3);
ax.YTick = 1:d(2);
grid(ax,'on');
ax.Box = 'on';
ax.YDir = 'reverse';
ax.YTickLabel = '';
ax.XTickLabel = '';
ax.FontSize = 7;
pl = find(pl);
for p=pl
    switch p
        case 1;ax.YTickLabel = labels{1,1};
        case 2;ax.XTickLabel = labels{1,2};
        case 3
            ax1 =  axes('position',[0.825,0.2,0.025,0.6]);
            imagesc(ax1,[0:100]');
            ax1.YDir = 'normal';
            colormap(ax1,gmap);
            ax1.XTick = [];
            ax1.YAxisLocation ='right';
    end
end

end
%------------- END OF CODE --------------
