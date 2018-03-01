function [im,a]=vane(w,gp,cmap,fig,thr,pw,ph,label)
%% PLOT.HIVE: One line description of what the function or script performs
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
%  VERSION:  0.0 CREATED: 19-Jul-2017 13:16:09
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
hive(input01,input02,input03,input04)
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
if ~exist('cmap','var');cmap = jet(nc);end
if ~exist('fig','var');fg= figure();else;fg = figure(fig);end
if ~exist('pw','var');pw = 500;end
if ~exist('ph','var');ph = 500;end
if ~exist('label','var');label = 1;end

set(fg,'units','pixels','position',[0 0 pw ph],'Visible','on','resize','off')


if ~exist('thr','var');thr= 0;end
nn=height(gp);
A = zeros(nn);
idx = find(triu(true(nn),1));
A(idx)=w;
A(A<thr)=0;
nc = numel(unique(gp.cid));
ang=0.5*pi:2*pi/nc:2.5*pi;
lang=0.47*pi:2*pi/nc:2.47*pi;
X=cos(ang);
Y=sin(ang);
x=cos(lang);
y=sin(lang);
labels = categories(gp.name);
nxy = zeros(nn,4);
clf
for ii=1:nc
    ix = find(gp.name==labels{ii});
    ld = [linspace(X(ii),X(ii)*.25,numel(ix))',linspace(Y(ii),Y(ii)*.25,numel(ix))'];
    line([X(ii)*.25,X(ii)],[Y(ii)*.25,Y(ii)]);
    for jj=1:numel(ix)
        nxy(ix(jj),:) = [ld(jj,1),ld(jj,2),ii,jj];
    end  
    if label
    t = atan2(Y(ii),X(ii));
    if abs(t) > pi/2        
        text(x(ii),y(ii),labels{ii},'rotation', 180*(t/pi + 1),'HorizontalAlignment','left','interpreter','none')
    else
        text(x(ii),y(ii),labels{ii},'rotation', t*180/pi,'HorizontalAlignment','right','interpreter','none')
    end
    end

end
t = linspace(0,1,100)';
f = @(t,xy1,cxy,xy2) kron((1-t).^2,xy1) + kron(2*(1-t).*t,cxy) +  kron(t.^2,xy2); % define a second order Bernstein polynomial
EW = c3nl.scale(A,1,4);
for ii=1:nn
    J = find(A(ii,:));
    if ~isempty(J)       
        x1 = nxy(ii,1);y1 = nxy(ii,2);
        for jj=1:numel(J)           
            x2 = nxy(J(jj),1);y2 = nxy(J(jj),2);
            ew = EW(ii,J(jj));
            if nxy(J(jj),3)==nxy(ii,3)
                [x,y]=gencurve(x1,y1,x2,y2,1);
                patch([real(x),NaN],[real(y),NaN],'k','linewidth',ew,'edgealpha',0.85,'edgecolor',cmap(nxy(J(jj),3),:)); % same cluster vector on the right
            else 
                
                xy=f(t,[x1,y1],[0,0],[x2,y2]);
                x = xy(:,1);y=xy(:,2);
%                x = linspace(x1,x2,100);y = linspace(y1,y2,100);
                patch([x;NaN],[y;NaN],'k','linewidth',ew*0.5,'edgealpha',0.5,'edgecolor',[0.5,0.5,0.5]); % same cluster vector on the right
            end
        end
    end
end
disp('')
deg = c3nl.scale(sum(apply.symAdj(A>0,'upper')),0.005,0.05);
for ii=1:size(nxy,1)
        [cx,cy]=plot.circle(nxy(ii,1),nxy(ii,2),deg(ii));
       try
        patch(cx,cy,'k','FaceAlpha',1,'facecolor',cmap(nxy(ii,3),:));% plot the nodes ,'parent', ax
       catch
          keyboard 
       end
end
axis image
axis 'off'
tmp = getframe(fg);
im = tmp.cdata;
a = rgb2gray(tmp.cdata)~=240;

end

function [x,y]=gencurve(x1,y1,x2,y2,rl)
pc = [max([x1,x2])-range([x1,x2])/2,max([y1,y2])-range([y1,y2])/2];
r =  sqrt((x2-pc(1))^2+(y2-pc(2))^2);
d = sqrt((x2-x1)^2+(y2-y1)^2); % Distance
ar = atan2(-(x2-x1),y2-y1);
al = atan2(x2-x1,-(y2-y1)); % Perpendicular bisector angle
b = asin(d/2/r); % Half arc angle
e = sqrt(r^2-d^2/4); % Distance, center to midpoint
cl = linspace(al-b,al+b); % Arc angle range
cr = linspace(ar-b,ar+b); % Arc angle range
xl = (x1+x2)/2-e*cos(al)+r*cos(cl); % Cartesian coords. of arc
yl = (y1+y2)/2-e*sin(al)+r*sin(cl);
xr = (x1+x2)/2-e*cos(ar)+r*cos(cr); % Cartesian coords. of arc
yr = (y1+y2)/2-e*sin(ar)+r*sin(cr);

if rl 
    if (x1 - x2)*(yl(50) - y2) > (y1 - y2)*(xl(50) - x2)% if this is true it means a flip
        x = xl;y = yl;
    else 
        x = xr;y = yr;
    end
end
    if (x1 - x2)*(yr(50) - y2) > (y1 - y2)*(xr(50) - x2)
        x = xr;y = yr;
    else 
        x = xl;y = yl;
    end
end




%------------- END OF CODE --------------
