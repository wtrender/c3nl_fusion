function h = projectedAdj(adj,T,obj,fig)
%% PLOT.PROJECTEDADJ: takes an undirected weighted graph and projects it onto the brain
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
%  VERSION:  0.0 CREATED: 19-Jun-2017 09:10:16
%
%% INPUTS:
%    adj - the undirected graph (size n x n)
%    T - the node infomration table (n long containing mni cordinates,
%    numeric id and labels)
%    obj - a structure containing the follwoing fields
%        or - the brain oreintation needed
%        thr - a limiting threshold to apply to the weights
%        side - 'positive' or 'negative'
%        cmap - the color map to use
%        bkg -
%        layers -
%        nodes -
%
%
%
%
%% OUTPUT:
%   h a figure handle
%
%% EXAMPLES:
%{
h = plot.projectedAdj(adj,T,obj)
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
if ~exist('fig','var');fig=1;end
output = setparams(adj,T,obj);
figure(fig);clf
for ii=1:5
    id = ['p' num2str(ii)];
    ax.(id) = axes('position',[0 0 1 1]);axis(ax.(id), 'off');
end

imagesc(output.Alpha,'AlphaData',output.Alpha,'parent',ax.p1);
colormap(ax.p1,output.bkg);axis(ax.p1, 'image');axis(ax.p1, 'off');
imagesc(output.im.img,'AlphaData',output.im.alpha,'parent',ax.p2);
axis(ax.p2, 'image');axis(ax.p2, 'off');
imagesc(output.AO,'AlphaData',(output.AOa.*output.Alpha*.65),'parent',ax.p3);
colormap(ax.p3,bone);axis(ax.p3, 'off');axis(ax.p3, 'image');
cla(ax.p4)
cla(ax.p5)
if output.nodes
    %labels = output.nl(output.nix);
    num = num2cell(output.nid);
    [yy,xx,zz]=c3nl.deal(output.xyz);
    v = output.cortex.V;
    switch output.or
        case 'a'
            [XY,ix]=unique([xx,yy],'rows');
            XY= flipud(XY);
            ix= flipud(ix);
            X = XY(:,1);Y = X;Z = zz(ix);
            Y(yy(ix)<0) = round(min(v(:,2)));
            Y(yy(ix)>=0) = round(max(v(:,2)));
            X = linspace(max(v(:,1)),min(v(:,1)),numel(X));
            for ii=1:numel(Y)
                line('XData',[xx(ix(ii));X(ii);X(ii)],'YData',[yy(ix(ii));yy(ix(ii));Y(ii)],'ZData',[min(v(:,3))*[1;1;1]],'color','k','linewidth',2,'LineStyle',':','parent',ax.p4);
                [cx,cy]=plot.circle(X(ii),Y(ii),6);
                patch('XData',cx,'YData',cy,'ZData',min(v(:,3))*cx.^0+10,'facecolor',[0.75,0.75,0.75],'edgecolor','none','parent',ax.p4);
                line('XData',cx,'YData',cy,'ZData',min(v(:,3))*cx.^0+10,'color','k','linewidth',1,'parent',ax.p4)
                text(X(ii),Y(ii),Z(ii),num(ix(ii)),'parent',ax.p5,'FontUnits', 'normalized','fontsize',0.04,'HorizontalAlignment','center','VerticalAlignment','middle','color','k','fontweight','bold');
            end
        case {'rs','ls'}
            [X,ix] = sort(xx,'descend');
            Y = yy;Z = zz;X = xx;
            Z(zz(ix)<0) = round(min(v(:,3)));
            Z(zz(ix)>=0) = round(max(v(:,3)));
            X = linspace(max(v(:,1)),min(v(:,1)),numel(Z));
            for ii=1:numel(yy)
                line('XData',[xx(ix(ii));X(ii);X(ii)],'YData',[min(v(:,3))*[1;1;1]],'ZData',[zz(ix(ii));zz(ix(ii));Z(ii)],'color','k','linewidth',2,'LineStyle',':','parent',ax.p4);
                [cx,cy]=plot.circle(X(ii),Z(ii),6);
                patch('XData',cx,'YData',min(v(:,3))*cx.^0+10,'ZData',cy,'facecolor',[0.75,0.75,0.75],'edgecolor','none','parent',ax.p4);
                line('XData',cx,'YData',min(v(:,3))*cx.^0+10,'ZData',cy,'color','k','linewidth',1,'parent',ax.p4)

                text(X(ii),Y(ii),Z(ii),num(ix(ii)),'parent',ax.p5,'FontUnits', 'normalized','fontsize',0.04,'HorizontalAlignment','center','VerticalAlignment','middle','color','k','fontweight','bold');
            end
    end
    daspect(ax.p4,[1 1 1]);axis(ax.p4, output.im.lim);daspect([1 1 1]);view(ax.p4,output.az,output.el);ax.p4.YDir='reverse';
    daspect(ax.p5,[1 1 1]);axis(ax.p5, output.im.lim);daspect([1 1 1]);view(ax.p5,output.az,output.el);ax.p5.YDir='reverse';
end

end


function output = setparams(adj,T,obj)

prop = {'or','side','bkg','thr','cmap','layers','el','az','nodes','node_labels','node_id','T','adj','cortex','lim','output'};
if ~isstruct(obj);ia = zeros(numel(prop),1); else; ia = ismember(prop,fieldnames(obj));end
for ii=1:length(ia)
    tmp = prop{ii};
    switch tmp
        case 'or';if ia(ii);output.or = obj.(tmp);else;output.or='a';end
        case 'side';if ~ia(ii);side='positive';end
        case 'nodes';if ia(ii);output.nodes = obj.(tmp);else; output.nodes =false;end
        case 'bkg';if ~ia(ii);bkg=[235,245,255]./256;end
        case 'thr';if ia(ii);thr = obj.(tmp);else;thr = 0;end
        case 'cmap'
            if ia(ii);cmap = obj.(tmp);
            else
                hmap = hsv(360);
                switch side
                    case 'positive'
                        cmap = [0.5,0.5,0.5;hmap(15,:).^.5;hmap(15,:)];
                        clim = round(max(adj(:))).*[0,1];
                    case 'negative'
                        cmap = [hmap(220,:);hmap(220,:).^.5;0.5,0.5,0.5];
                        clim = round(min(adj(:))).*[1,0];
                end
            end
        case 'layers';if ia(ii);layers = obj.(tmp);else; load atlas/layers96.mat;end
        case 'az';if ~ia(ii);az = layers.(output.or).az;end
        case 'el';if ~ia(ii);el = layers.(output.or).el;end
        case 'lim';    lim=[-120,80,-80,80,-80,100];
        case 'T'
            n = height(T);
            ix=find(triu(true(n),1));
        case 'adj'
            tmp = zeros(n,n,size(adj,2));
            for m=1:size(adj,2)
                tmp1 = zeros(n);
                tmp1(ix(adj(:,m)~=0)) = adj(adj(:,m)~=0,m);
                tmp(:,:,m) = apply.symAdj(tmp1,'upper');
            end
            adj = tmp;
            switch side
                case 'positive';adj(adj<thr)=0;
                case 'negative';adj(adj>thr)=0;
            end
        case 'node_labels';if ia(ii);node_labels = obj.(tmp);else;node_labels = matlab.lang.makeUniqueStrings(T.name_AAL2);end
        case 'node_id';if ia(ii);node_id = obj.(tmp);else;node_id = T.ROIid;end
        case 'cortex';if ia(ii);layers = obj.(tmp);else; load atlas/cortex.mat;end
        case 'output'
            output.Alpha = layers.(output.or).Alpha;
            output.AO = layers.(output.or).ao.img;
            output.AOa = layers.(output.or).ao_alpha.img;
            ph=fix((108*20)/2.54);
            [output.im,nodes] = genImage(adj,cortex,T,az,el,cmap,lim,ph);
            output.bkg = bkg;
            output.nl = nodes.AAL;
            output.nid = nodes.id;
            output.xyz = nodes.xyz;
            output.az =az;
            output.el =el;
            output.cortex = cortex;
    end
end
end

function [output,nodes] = genImage(adj,cortex,T,az,el,cmap,lim,ph)

nw = squeeze(sum(adj));
if numel(size(adj))==2
    nw = nw(:);
end
if nnz(nw~=0)
    nodes = table();
    clim = round(max(abs(nw(:)))).*[0,1];
    for idx=find(nw)'
        [I,J]= ind2sub(size(nw),idx);
        tmap = [cmap(J,:).^.5;cmap(J,:);cmap(J,:).^2];
        nc = c3nl.ind2Cmap(tmap,[clim(:);nw(idx)]);
        try
            tmp = table(I,T.name_AAL2(I),[T.X(I,:),T.Y(I,:),T.Z(I,:)],nc(end,:),nw(idx),J,'VariableNames',{'id','AAL','xyz','cmap','nw','class'});
        catch 
            keyboard
        end
        nodes = [nodes;tmp];
    end
    links = table();
    A = sum(adj>0,3);
    A(tril(true(size(A)),-1)) = 0;
    [ixI,ixJ] = find(A);
    llim = round(max(adj(:))).*[0,1];
    for ii=1:numel(ixJ)
        J = find(adj(ixI(ii),ixJ(ii),:));
        for m=1:numel(J)
            tmap = [cmap(J(m),:).^.5;cmap(J(m),:);cmap(J(m),:).^2];
            lc = c3nl.ind2Cmap(tmap,[llim(:);adj(ixI(ii),ixJ(ii),J(m))]);
            Centroid1= [T.X(ixI(ii),:),T.Y(ixI(ii),:),T.Z(ixI(ii),:)];
            Centroid2= [T.X(ixJ(ii),:),T.Y(ixJ(ii),:),T.Z(ixJ(ii),:)];
            tmp = table(ixI(ii),ixJ(ii),Centroid1,Centroid2,lc(end,:),adj(ixI(ii),ixJ(ii),J(m)),J(m),'VariableNames',{'n1','n2','xyz1','xyz2','cmap','w','class'});
            links = [links;tmp];
        end
    end
    
    fg = figure(99);clf;
    set(fg,'units','pixels','position',[0 0 ph ph],'Visible','on','resize','off')
    ax=axes('position' ,[0 0 1 1]);axis off;
    patch('Vertices',cortex.V,'Faces',cortex.F,'FaceColor','r','FaceAlpha',0,'EdgeColor','none');% create the axis based on the cortex
    hold on
    w = c3nl.scale(links.w,.5,2.5);
    for ii=1:numel(w) % plot the links
        tmp = [links.xyz1(ii,[2,1,3]);links.xyz2(ii,[2,1,3])];
        [X,Y,Z] = get.tube(tmp',15,w(ii));
        surf(ax,X,Y,Z,'FaceLighting','gouraud','facecolor',links.cmap(ii,:),'edgecolor','none')
    end
    w = c3nl.scale(nodes.nw,2.5,6.5);
    for ii=1:numel(w) % go over the nodes
       [X,Y,Z] = ellipsoid(nodes.xyz(ii,2),nodes.xyz(ii,1),nodes.xyz(ii,3),w(ii),w(ii),w(ii),20);           
       surf(ax,X,Y,Z,'FaceLighting','gouraud','facecolor',nodes.cmap(ii,:),'faceAlpha',0.75,'edgecolor','none'); 
    end
    view(az,el);axis 'off';daspect([1 1 1]);
    camlight(-80,-10); camlight(80,-10);
    ax.YDir = 'reverse';  axis(lim)  
    tmp = getframe(gcf);
    output.img = tmp.cdata;
    output.lim = [xlim,ylim,zlim];
    output.alpha = rgb2gray(tmp.cdata)~=240;
end
end

%------------- END OF CODE --------------
