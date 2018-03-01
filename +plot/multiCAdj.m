function multiCAdj(T,W,cmap,fs,fig,side)
%% PLOT.MULTICADJ: One line description of what the function or script performs
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
%  VERSION:  0.0 CREATED: 13-Jun-2017 14:38:29
%
%% INPUTS:
%    nodes - node labels
%    W - a edge weight vector/matrix of (nodes*nodes-1)/2 x classes
%    mode - 'positive', 'negative','both' the sign of the wieghts to plot
%    cmap - a color map with [R G B] x classes
%    fig - figure number
%    fs  - font size
%
%
%% OUTPUT:
%   weight -
%   degree -
%   A -
%% EXAMPLES:
%{
[weight,degree,A]= plot.multiCAdj(nodes,W,mode,cmap)
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

nn = height(T);
if exist('fig','var')&&~isempty(fig);fg = figure(fig);else;fg = figure;end
if ~exist('thr','var');thr = 0;end
if ~exist('cmap','var');cmap = jet;end
if ~exist('fs','var');fs = 12;end
if ~exist('side','var');side= 'positive';end

%% convert fs upper triangle weights to sparse form

ix=find(triu(true(nn),1));
tmp = zeros(nn,nn,size(W,2));
for m=1:size(W,2)
    tmp1 = zeros(nn);
    tmp1(ix(W(:,m)~=0)) = W(W(:,m)~=0,m);
    tmp(:,:,m) = apply.symAdj(tmp1,'upper');
end
adj = tmp;
switch side
    case 'positive';adj(adj<0)=0;lim =[0,1];
    case 'negative';adj(adj>0)=0;lim =[0,1];
    otherwise;lim =[-1,1];
end

nw = squeeze(sum(adj~=0)); % degree
if size(adj,3)==1;nw=nw';end
if nnz(nw~=0)
    n = nnz(any(nw,2));
    rad = linspace(0.5*pi,2.5*pi,n+1);
    X=cos(rad);
    Y=sin(rad);
    t = linspace(0,1,101);
    nodes = table();
    clim = round(max(abs(nw(:)))).*[0,1];
    c = 0;
    nid = zeros(n,1);
    for ii=1:nn
        ix = find(nw(ii,:));
        if ~isempty(ix)
            c = c+1;
            nid(c)=ii;
            for J=ix
                tmap = [cmap(J,:).^.5;cmap(J,:);cmap(J,:).^2];
                nc = c3nl.ind2Cmap(tmap,[clim(:);nw(ii,J)]);
                tmp = table(ii,T.name_AAL2(ii),X(c),Y(c),nc(end,:),nw(ii,J),J,'VariableNames',{'id','AAL','X','Y','cmap','nw','class'});
                nodes = [nodes;tmp];
            end
        end
    end
    links = table();
    A = sum(adj~=0,3);
    A(tril(true(size(A)),-1)) = 0;
    [ixI,ixJ] = find(A);
    llim = round(max(adj(:))).*lim;
    for ii=1:numel(ixJ)
        J = find(adj(ixI(ii),ixJ(ii),:));
        for m=1:numel(J)
            tmap = [cmap(J(m),:).^.5;cmap(J(m),:);cmap(J(m),:).^2];
            lc = c3nl.ind2Cmap(tmap,[llim(:);abs(adj(ixI(ii),ixJ(ii),J(m)))]);
            ix1 = find(nid==ixI(ii));ix2 = find(nid==ixJ(ii));
            tmp = table(ixI(ii),ixJ(ii),rad(ix1),rad(ix2),lc(end,:),adj(ixI(ii),ixJ(ii),J(m)),J(m),'VariableNames',{'n1','n2','rad1','rad2','cmap','w','class'});
            links = [links;tmp];
        end
    end
else 
    return
end

% plot links
clf
ax = axes('position',[0.25 0.25 0.5 0.5]);
hold on
clim = round(max(abs(links.w(:)))).*[0,1];
sw = c3nl.scale([abs(links.w);clim(:)],1,5);
sa = c3nl.scale([abs(links.w);clim(:)],0.5,0.85);

for jj=1:height(links)
    [x,y] = plot.arc(links.rad1(jj),links.rad2(jj));
    patch([x(:);NaN],[y(:);NaN],'k','linewidth',sw(jj),'edgealpha',sa(jj),'edgecolor',links.cmap(jj,:),'parent', ax);
end
% Plot the unit circle
[cx,cy]=plot.circle(0,0,1);% get a unit circle
patch([cx';NaN],[cy';NaN],'k','linewidth',3,'edgealpha',0.5,'edgecolor','k','parent', ax);

clim = round(max(abs(nodes.nw(:)))).*[0,1];
sw = c3nl.scale([nodes.nw;clim(:)],0.02,0.06);
labels = arrayfun(@(a,b) sprintf('%i. %s',a,b{1}), nodes.id,nodes.AAL,'un',0);
for ii=1:height(nodes)
    yy = nodes.Y(ii);
    xx = nodes.X(ii);
   [cx,cy]=plot.circle(xx,yy,sw(ii));% get small circles
    patch(cx,cy,'k','FaceAlpha',1,'facecolor',nodes.cmap(ii,:),'parent', ax);% plot the nodes 
    t = atan2(yy,xx);
    if abs(t) > pi/2
        text(1.1*xx,1.1*yy,labels{ii},'rotation', 180*(t/pi + 1),'fontsize',fs,'HorizontalAlignment','right','interpreter','none','parent', ax)
    else
        text(1.1*xx,1.1*yy,labels{ii},'rotation', t*180/pi,'fontsize',fs,'HorizontalAlignment','left','interpreter','none','parent', ax)
    end
end

axis image
axis([1.2*xlim,1.2*ylim])

axis 'off'

end
%------------- END OF CODE --------------
