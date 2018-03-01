function [measuresOrder,featuresOrder,featureClusters,measuresClusters,fhsv]=HCLvert(X,obj,dmethod,lmethod,cmap,mmap)

disp('')
set(0,'RecursionLimit',length(X));
[m,n] = size(X);
fs = 14;
[Title,xlabel,Y1,Y2,Y3,Y4,Yc1,Yc2,Yc3,Yc4,c0,c1,c2,c3,clim,fg,ax,zs,srt] = getParams(X,obj);
if zs ; Data = [zscore(X)]'; else Data = X';end
if ~exist('dmethod','var');dmethod='euc';end
if ~exist('lmethod','var');lmethod='weighted';end
Z_features = linkage(pdist(Data,dmethod),lmethod);%weighted,average,complete ,'centroid','median','single','ward'
Z_measures = linkage(pdist(Data',dmethod),lmethod);
Y = inconsistent(Z_features,3);
[~,id]=sort(Y(:,4),'descend');
%[~,id]=max(Y(:,2))
featuresOrder = optimalleaforder(Z_features,pdist(Data));
axes(ax.dendoleft);[h,~,Z_features_order,lcf] = plot.dendrogram(Z_features,0,'orient', 'left', 'colorthreshold',Z_features(id(3),3),'Reorder',featuresOrder);
set(h, 'linewidth', 1);axis(ax.dendoleft,'off')
conn = (Z_features(:,3) < Z_features(id(3),3));

[featureClusters,L] = get.labeltree(Z_features, conn,lcf);
Y = inconsistent(Z_measures,1);
[~,id]=sort(Y(:,4),'descend');
measuresOrder = optimalleaforder(Z_measures,pdist(Data'));

% [lc{1:numel(h)}] = h.Color;
% lc = cell2mat(lc');
% lcl = double(categorical(cellstr(char(lc*10+57))))-1;
if ~exist('cmap','var')
    fhsv = hsv(max(L)-1);
else 
    fhsv = cmap;
end
%[lcf,featureClusters] = branch2node(lcf,featureClusters);
for ii=1:numel(lcf)
    if lcf(ii)
        h(ii).Color =  fhsv(lcf(ii),:);
        h(ii).LineWidth = 2;
    end
end
%0.7*max(Z_measures(:,3)) Z_measures(id(2),3)
axes(ax.dendotop);[h,~,Z_measures_order,lc]=dendrogram(Z_measures,0,'Reorder',measuresOrder, 'colorthreshold',Z_measures(end-1,3));
conn = (Z_measures(:,3) < Z_measures(id(1),3));

[measuresClusters,L] = get.labeltree(Z_measures, conn,lc);
%set(h, 'Color', 'k');
set(h, 'linewidth', 1);axis(ax.dendotop,'off')
imagesc(flipud(Data(featuresOrder,measuresOrder)),'parent',ax.heatmap);colormap( ax.heatmap,c0);axis(ax.heatmap,'off');
cb = linspace(clim(1),clim(2),11);
imagesc((cb)','parent',ax.cb)
colormap( ax.cb,c0);
set(ax.heatmap,'CLim',clim);
set(ax.cb,'CLim',clim);
ax.cb.YAxisLocation = 'right';
ax.cb.YTick = 1:11;
ax.cb.YTickLabel = num2cell((round(cb(ax.cb.YTick),2)));
ax.cb.XTick ='';
text(ax.cb,0.5,0,'$$\sum\beta$$','Interpreter','latex')
ax.cb.YDir = 'normal';
ax.cb.FontSize = fs;
set(ax.heatmap, 'Xticklabel', [], 'yticklabel',[]);
if ~isempty(xlabel)
text(ax.labelsright, repmat(0.01,n,1), ax.dendoleft.YTick,xlabel(Z_features_order),...
    'horizontalalignment','left','FontSize',fs,'interpreter','none')
set(ax.labelsright, 'ylim', [0.5,length(Z_features_order)+0.5], 'Visible', 'off');
end
gp = [~isempty(Y1),~isempty(Y2),~isempty(Y3),~isempty(Y4)];
gc = get.cmap([120 34 102;230 210 230]);

for ii=1:sum(gp)
    axl = ax.(['l' num2str(ii)]);axb = ax.(['b' num2str(ii)]);
    switch ii
        case 1;gt = Y1;str=[char(65) '. ' Yc1.type];cat = Yc1.cat;
        case 2;gt = Y2;str=[char(66) '. ' Yc2.type];cat = Yc2.cat;
        case 3;gt = Y3;str=[char(67) '. ' Yc3.type];cat = Yc3.cat;
        case 4;gt = Y4;str=[char(68) '. ' Yc4.type];cat = Yc4.cat;
    
    end
    imagesc(flipud(gt(measuresOrder)'),'parent',axb);
    text(axb, -round(m*0.1),1 ,str,...
        'horizontalalignment','left','FontSize',fs,'interpreter','none','fontweight','bold');
    axis(axb,'off');
    if ii>1;colormap(axb,gc);else; colormap(axb,gc);end
if ~isempty(cat)
    GC = c3nl.ind2Cmap(gc,double(cat));
    if numel(cat)<4
        [x,y] = get.fancyGrid(numel(cat),0.1,0.1,'same');
        y(:) = 0.5;
    else
        r = round(numel(cat)^.5);
        A = ones(r,1)*(numel(cat)/r);
        [x,y] = get.fancyGrid(A,0.1,0.1,'same');
    end
    x=flipud(x);
    for jj=1:numel(cat)
        text(axl,x(jj),y(jj),char(9608),'color',GC(jj,:),'horizontalalignment','left','FontSize',fs)
        text(axl,x(jj),y(jj),sprintf('    %s',char(cat(jj))),'FontSize',fs,'horizontalalignment','left','verticalalignment','middle');
    end
end
    axis(axl,[0,1,0,1]);axis(axl,'off');
end

set(ax.dendoleft, 'ylim', [0.5,length(featuresOrder)+0.5], 'Visible', 'Off');
set(ax.dendotop, 'xlim', [0.5,length(measuresOrder)+0.5], 'Visible', 'Off');

end

function [lc,featureClusters] = branch2node(lc,groups)
    % nodes
    [a1,b1]= hist(featureClusters(featureClusters~=0),unique(featureClusters(featureClusters~=0)));
    [~,ix1]= sort(a1); % sort by number of nodes in each group
    % brances
    [a2,b2]= hist(lc(lc~=0),unique(lc(lc~=0)));
    [~,ix2]= sort(a2); % sort by number of brances in each group
    tmplc = zeros(size(lc));
    tmpfeatureClusters = zeros(size(featureClusters));
    if numel(ix1)~=numel(ix2);disp('this function requires equal amount of groups across both nodes and brnaces');
    else
       for ii=1:numel(ix1)
          % go over both lists and make the id's match
          tmpfeatureClusters(featureClusters==b1(ix1(ii)))=ii;
          tmplc(lc==b2(ix2(ii)))=ii;
           
       end
    end
    lc = tmplc;
    featureClusters = tmpfeatureClusters;
end


function [Title,xlabel,Y1,Y2,Y3,Y4,Yc1,Yc2,Yc3,Yc4,c0,c1,c2,c3,clim,fg,ax,zs,srt] = getParams(X,obj)

prop = {'Title','xlabel','Y1','Y2','Y3','Y4','Yc1','Yc2','Yc3','Yc4','c0','c1','c2','c3','clim','fn','fg','ax','srt','zs'};
if ~isstruct(obj);ia = zeros(numel(prop),1); else [ia,~] = ismember(prop,fieldnames(obj));end
for ii=1:length(ia)
    tmp = prop{ii};
    switch tmp
        case 'Title';if ia(ii);Title=obj.(tmp);else Title=[];end
        case 'xlabel';if ia(ii);xlabel=obj.(tmp);else xlabel=[];end
        case 'Y1';if ia(ii);Y1=obj.(tmp);else Y1=[];end
        case 'Y2';if ia(ii);Y2=obj.(tmp);else Y2=[];end
        case 'Y3';if ia(ii);Y3=obj.(tmp);else Y3=[];end
        case 'Y4';if ia(ii);Y4=obj.(tmp);else Y4=[];end
        case 'Yc1';if ia(ii);Yc1=obj.(tmp);else Yc1=[];end
        case 'Yc2';if ia(ii);Yc2=obj.(tmp);else Yc2=[];end
        case 'Yc3';if ia(ii);Yc3=obj.(tmp);else Yc3=[];end
        case 'Yc4';if ia(ii);Yc4=obj.(tmp);else Yc4=[];end
        case 'c0';if ia(ii);c0=obj.(tmp);else c0=get.cmap(flipud([255 205 155;255 180 104;255 155 54;255 130 4;191 98 3;51 51 51;14 83 164;49 92 143;62 117 182;110 152 200;207 221 237]./255),101);end
        case 'c1';if ia(ii);c1=obj.(tmp);else c1=[];end
        case 'c2';if ia(ii);c2=obj.(tmp);else c2=[];end
        case 'c3';if ia(ii);c3=obj.(tmp);else c3=[];end
        case 'zs';if ia(ii);zs=obj.(tmp);else zs=0;end
        case 'srt';if ia(ii);srt=obj.(tmp);else srt=1;end
        case 'clim';if ia(ii);clim=obj.(tmp);else clim= round(max(abs(X(:)))).*[-1,1];end
        case 'fn';if ia(ii);fn=obj.(tmp);else fn=3;end
        case 'fg';if ia(ii);fg=obj.(tmp);else fg=figure(fn);clf;end
        case 'ax'
            if ia(ii);ax=obj.(tmp);
            else
                clf;b1 = .15;l1=.075;w1 = .72;h1 =.68;
                ax.heatmap = axes('Position', [l1 b1 w1 h1]);%[left bottom width height]
                ax.dendotop = axes('Position',  [l1 b1+h1 w1 0.16]);%axis('off');
                ax.dendoleft = axes('Position',  [0.01 b1 l1-.01 h1]);%axis('off');
                if ~isempty(xlabel); ax.labelsright = axes('Position',  [l1+w1 b1 0.15 h1]);end%axis('off');
                ax.cb = axes('Position',  [w1+l1+.15 b1 0.01 h1]);%axis('off');
                gp = [~isempty(Y1),~isempty(Y2),~isempty(Y3),~isempty(Y4)];
                if any(gp)
                    A = ones(sum(gp),1);
                    [x,y,w,h]=get.fancyGrid(A,0.05,0.001,'same');
                    
                    for jj=1:numel(x)
                        ax.(['b' num2str(jj)]) = axes('Position',  [l1 0.01+y(jj)*(b1-.05) w1 h*(b1-0.05)]);
                        ax.(['l' num2str(jj)]) = axes('Position',  [l1+w1 0.01+y(jj)*(b1-.05) w*.2 h*(b1-0.05)]);
                    end
                end
            end
            
    end
end
end


