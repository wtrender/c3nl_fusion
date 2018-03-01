function [excmat,ACC,F1,G,Prc,Rec,fg] = confusiongrad(Y, Yfit, class, r_title ,w, h,fc,fig,rmlabels)
%% PLOT.CONFUSIONGRAD: One line description of what the function or script performs
% [excmat,ACC,F1,G,Prc,Rec,fg] = confusiongrad(Y, Yfit, class, r_title ,w, h,fc,fig,rmlabels)
%   __           _
%  / _|         (_)
% | |_ _   _ ___ _  ___  _ __
% |  _| | | / __| |/ _ \| `_ \    :- Functional and Structural
% | | | |_| \__ \ | (_) | | | |      Integration of Neuroimages
% |_|  \__,_|___/_|\___/|_| |_|
%
%
% AUTHOR:  Eyal Soreq
%  EMAIL:  e.soreq14@imperial.ac.uk
%  AFFILIATION:  Imperial College London
%  VERSION:  0.0 CREATED: 10-Jun-2017 09:06:20
%
% INPUTS:
%    input01 -
%    input02 -
%    input03 -
%    input04 -
%
%
% OUTPUT:
%
% EXAMPLES:
%{
confusiongrad(input01,input02,input03,input04)
%}
%
% DEPENDENCIES:
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
if ~exist('class','var')||isempty(class);class=categories(Yfit);end
if ~exist('r_title','var')||isempty(r_title);r_title='Confusion Matrix';end
if ~exist('w','var')||isempty(w);w=600;end
if ~exist('h','var')||isempty(h);h=600;end
if ~exist('fc','var')||isempty(fc);fc=0.05;end
if ~exist('fig','var')||isempty(fig);fig=1;end
if ~exist('rmlabels','var')||isempty(rmlabels);rmlabels = zeros(4,1);end
if ~exist('grad','var')||isempty(grad);grad = 0;end



if ~(iscategorical(Y)&&iscategorical(Yfit))
    elem = [numel(unique(Y)),numel(unique(Yfit))];
    if diff(elem)
        [~,id]=max(elem);
        switch id;case 1;vset = unique(Y);case 2;vset = unique(Yfit);end
    else ; vset = unique(Y);
    end
    try
        Y = categorical(Y,vset,class);
        Yfit = categorical(Yfit,vset,class);
    catch
        keyboard
    end
end
d = numel(class);
tp = eye(d)>0;
cmat = confusionmat(Yfit,Y,'ORDER',class)';
ACC = sum(Y == Yfit)/length(Y);ACC = [ACC, 1-ACC]*100;
% True rate the proportion between the true positive and the ground truth
% aka Sensitivity Recal
TR =  diag(cmat)./ sum(cmat,2);TR = [TR,1-TR]*100;
% Predictive Value the proportion between the true positive and the any
% value classified in class k aka Precision
PV = diag(cmat)./sum(cmat)';PV = [PV,1-PV]*100;
F1score =   2*((PV(:,1).*TR(:,1))./(PV(:,1) + TR(:,1)));
Prc = nanmean(PV(:,1));
Rec = nanmean(TR(:,1));
% can be also calculated as F1score = 2./(1./PV(:,1)+1./TR(:,1));
Gscore = (PV(:,1).*TR(:,1)).^0.5;
F1 = nanmean(F1score);
G = nanmean(Gscore);

if exist('r_title','var')
    fg = figure(fig);
    set(fg,'units', 'pixels' ,'pos', [100 100 w h]);
    %colormap([115 231 179;232 130 130;128 128 128; 119 148 233]./255)
    colormap([230 230 230;180 180 180;140 140 140; 119 148 233]./255)
    clf
    m = numel(class);
    us = 0.65/(m+1);
    
    ax1 = axes('Position', [0.2, 0.15+us, us*m, us*m]);% [left bottom width height]
    ax2 = axes('Position', [0.2+us*m, 0.15+us, us, us*m]);% [left bottom width height]
    ax3 = axes('Position', [0.2, 0.15, us*m, us]);% [left bottom width height]
    ax4 = axes('Position', [0.2+us*m, 0.15, us, us]);% [left bottom width height]
    
    cmap1=get.cmap([0.8 0.1 0.1;0.6,0.2,0.2;0.8 0.8 0.8;0.2 0.6 0.2;0 1 0],101);
    cmap2=get.cmap([1,1,1;0.8 0.8 0.8;0.2 0.6 0.2;0 1 0],101);
    excmat = round(100*cmat/numel(Y),2);
    excmat(~eye(numel(class))) = -excmat(~eye(numel(class)));
    imagesc(excmat,'parent',ax1)
    colormap(ax1,cmap1)
    set(ax1,'CLim',[-100/numel(class),100/numel(class)]);
    set(ax1,'xtick',.5:numel(class), 'ytick',.5:numel(class),...
        'Xticklabel','','Yticklabel','','GridLineStyle','-');
    grid on;
    imagesc(TR(:,1),'parent',ax2)
    colormap(ax2,cmap1)
    set(ax2,'CLim',[0,100]);
    set(ax2,'xtick','', 'ytick',.5:numel(class),...
        'Xticklabel','','Yticklabel','','GridLineStyle','-');
    grid on;
    imagesc(PV(:,1)','parent',ax3)
    colormap(ax3,cmap1)
    set(ax3,'CLim',[0,100]);
    set(ax3,'xtick',.5:numel(class), 'ytick','',...
        'Xticklabel','','Yticklabel','','GridLineStyle','-');
    grid on;
    imagesc(ax4,ACC(:,1)')
    colormap(ax4,cmap1)
    set(ax4,'CLim',[0,100]);
    set(ax4,'xtick',.5:1.5, 'ytick',.5:1.5,...
        'Xticklabel','','Yticklabel','','GridLineStyle','-');
    grid on;
    
    [p,q] = meshgrid(1:d+1, 1:d+1);
    
    
    excmat = zeros([size(cmat)+1,2]);
    excmat(1:d,1:d,1) = cmat;
    excmat(1:d,1:d,2) = round(100*cmat/numel(Y),2);
    excmat(end,1:d,1) = PV(:,1)';
    excmat(end,1:d,2) = PV(:,2)';
    excmat(1:d,end,1) = TR(:,1)';
    excmat(1:d,end,2) = TR(:,2)';
    excmat(d+1,d+1,:) = ACC;
    pairs = [p(:) q(:)];
    
    ax5 = axes('Position', [0.2, 0.15, us*(m+1), us*(m+1)]);
    if ~isempty(fc)
        for ii=1:length(pairs)
            if pairs(ii,1)<=d && pairs(ii,2)<=d
                text(pairs(ii,2),pairs(ii,1)-0.1,...
                    sprintf('%i',excmat(pairs(ii,1),pairs(ii,2),1)),...
                    'HorizontalAlignment','center','FontUnits', 'normalized','fontsize',fc)
                text(pairs(ii,2),pairs(ii,1)+0.1,...
                    sprintf('%3.2f%%',excmat(pairs(ii,1),pairs(ii,2),2)),...
                    'HorizontalAlignment','center','FontUnits', 'normalized','fontsize',fc)
            else
                text(pairs(ii,2),pairs(ii,1)-0.1,...
                    sprintf('%3.2f%%',squeeze(excmat(pairs(ii,1),pairs(ii,2),1))),...
                    'HorizontalAlignment','center','FontUnits', 'normalized','fontsize',fc,...
                    'color',[0 0.4 0])
                text(pairs(ii,2),pairs(ii,1)+0.1,...
                    sprintf('%3.2f%%',squeeze(excmat(pairs(ii,1),pairs(ii,2),2))),...
                    'HorizontalAlignment','center','FontUnits', 'normalized','fontsize',fc,...
                    'color',[0.7 0 0])
            end
        end
        axis([0.5 d+1.5 0.5 d+1.5])
        axis off
        ax5.YDir = 'reverse';
        if ~mod(d,2);sh=0.5;else sh=0;end
        ix = find(~rmlabels);% T,B L R
        for ii=1:numel(ix)
            switch ix(ii)
                case 1
                    text(ceil(d/2)+sh,0,'Predicted condition','HorizontalAlignment','center','FontUnits', 'normalized','fontsize',fc)
                    text(1:d,zeros(d,1)+.3,class,'HorizontalAlignment','left','FontUnits', 'normalized','fontsize',fc,'rotation',30)
                case 2
                    text(ceil(d/2)+sh,d+1.7,'Precision','HorizontalAlignment','center','FontUnits', 'normalized','fontsize',fc)
                case 3
                    text(d+1.7,ceil(d/2)+sh,'Sensitivity','Rotation',-90,'HorizontalAlignment','center','FontUnits', 'normalized','fontsize',fc)
                case 4
                    text(0,ceil(d/2)+sh,'Ground truth','Rotation',90,'HorizontalAlignment','center','FontUnits', 'normalized','fontsize',fc)
                    text(zeros(d,1)+.3,1:d,class,'Rotation',30,'HorizontalAlignment','right','FontUnits', 'normalized','fontsize',fc)
            end
        end
        text(0.5,d+2,r_title,'HorizontalAlignment','left','FontUnits', 'normalized','fontsize',fc)
        text(d+1,.3,sprintf('F1=%3.2f%%',F1),'HorizontalAlignment','center','FontUnits', 'normalized','fontsize',fc,'color',[.9,0,0],'fontweight','bold')
    end
else
    excmat = cmat;
    
end


end
%------------- END OF CODE --------------
