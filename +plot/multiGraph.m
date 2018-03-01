function multiGraph(T,W,gp,X,factors,cmap,fig,w,h)
%% PLOT.MULTIGRAPH: One line description of what the function or script performs
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
%  VERSION:  0.0 CREATED: 10-Jan-2018 04:20:20
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
multiGraph(input01,input02,input03,input04)
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
fs = any(W~=0,2);
if exist('fig','var')&&~isempty(fig);fg = figure(fig);else;fg = figure;end
if ~exist('w','var');w = 1200;end
if ~exist('h','var');h = 320;end
set(fg,'units','pixels','position',[0 0 w h],'Visible','on','resize','off');

if ~exist('cmap','var');cmap = jet;end
ix=find(triu(true(nn),1));
tmp = zeros(nn,nn,size(W,2));
for m=1:size(W,2)
    tmp1 = zeros(nn);
    tmp1(ix(W(:,m)~=0)) = W(W(:,m)~=0,m);
    tmp(:,:,m) = apply.symAdj(tmp1,'upper');
end
adj = tmp;
fc = categorical(factors{1});
ng = squeeze(any(adj~=0))*[1:2:numel(fc)*2-1]';
pc = [[1:numel(fc);1:numel(fc)]';nchoosek(1:numel(fc),2)];
pc = fc(pc(:,1)).*fc(pc(:,2));


nw = squeeze(sum(adj~=0)); % degree
if size(adj,3)==1;nw=nw';end
if nnz(nw~=0)
    n = nnz(any(nw,2));
    Aa = (sum(adj~=0,3));
    ix = any(Aa,2);
    NL = T.ROIid(ix);
    G = graph(Aa(ix,ix));
    g = plot(G,'Layout','force');
    NXY = [g.XData',g.YData']; % node xy 
    %NW(ix,fs)
    NG = categorical(ng(ix),[1,3,5,4,6,8],cellstr(pc));
    ND = sum(squeeze((sum(adj~=0))),2); % degree
    EP = G.Edges.EndNodes;
    EG = NG(EP(:,1));
    for ii=1:size(EP,1)
        if EG(ii)~=NG(EP(ii,2))
            p = cellfun(@(x) nnz(c3nl.strDetect([strsplit(char(NG(EP(ii,1)))),strsplit(char(NG(EP(ii,2))))],x)),cellstr(fc),'un',0);
            [~,id]=find(([p{:}]));
            EG(ii) = fc(id(1)).*fc(id(2));
        end
    end
    c = 0;
    LIM = [0.05,max(X(:))*0.8];
    [x,y,h,w]=get.fancyGrid(ones(numel(factors{1}),1),0.05,0.05,'same');
    clf;
    %cmap = ones(6,3).*0.3;
    im = cell(3,1);
    for jj=3:-1:1
        fg = figure(1);
        c= 0;
        %clf;
    for ii=jj:3:numel(gp)
        c= c+1;
        ew = X(ii,fs);
        ax = axes('position',[y(c),x(c),w,h]);
        name = [factors{1}{c}];
        plot.graph(NXY,ND(ix),NG,arrayfun(@(x) num2str(x),NL,'un',0),EP,ew,EG,LIM,repmat(cmap(jj,:),6,1),name,ax);
    end
    %tmp = getframe(fg);
    %im{jj} = rgb2gray(tmp.cdata);
    end
    %clf ;imagesc(c3nl.scale(double(cat(3,im{3},im{2},im{1}))));
    axis off;

end

%factors = [{unique(XX.domain)},{unique(XX.load)}];



%{
use graph to extract network layout 
transverse over the multiclass network to plot the edges as lines if only
one and as curves for multi-edge connection 
store each curve in an array 
duplicate the network across axes 
adapt node size and color based on degree
adapt edge color based on hirrchical link weight

%}

%------------- END OF CODE --------------
