function graph(nxy,nd,ng,nl,ep,ew,eg,lim,cmap,name,ax)

%% PLOT.GRAPH: One line description of what the function or script performs
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
%  VERSION:  0.0 CREATED: 11-Jan-2018 09:54:41
%
%% INPUTS:
%    nxy - node xy position
%    nd - node degree
%    ng - node group
%    nl - node label
%    ep - edge parents 
%    ew - edge weight 
%    eg - edge group 
%    cmap - color map per group 
%    ax - axis handle 
%
%% OUTPUT:
%
%% EXAMPLES:
%{
plot.graph(nxy,nd,ng,nl,ep,ew,eg,cmap,name,ax)
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
r = [(range(nxy(:))./numel(ng)).^1.6,(range(nxy(:))./numel(ng)).^1.3];
cla(ax);
eW = c3nl.scale([lim,ew],0.25,8);
eA = c3nl.scale([lim,ew],0.2,0.5);
eW(1:2)=[];
eA(1:2)=[];
for ii=1:size(ep,1)
   patch([nxy(ep(ii,1),1);nxy(ep(ii,2),1);NaN],[nxy(ep(ii,1),2);nxy(ep(ii,2),2);NaN],'k','linewidth',eW(ii),'edgealpha',1,'edgecolor',cmap(double(eg(ii)),:),'parent', ax);
end
nr = c3nl.scale(nd,r(1)*1.25,r(1)*1.25);
for ii=1:size(nxy,1)
   [xp,yp] = plot.circle(nxy(ii,1),nxy(ii,2),nr(ii));
   patch('XData',xp,'YData',yp,'FaceColor','w','parent',ax); 
   text(nxy(ii,1),nxy(ii,2),nl{ii},'Color','k','FontUnits','normalized','FontSize',0.035,'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','middle')
end
axis tight
axis off
x = xlim;
y = ylim;
text(x(1), y(1),name,'FontUnits','normalized','FontSize',0.05,'FontWeight','bold');


%------------- END OF CODE --------------
