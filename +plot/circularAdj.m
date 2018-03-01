function circularAdj(fig,T,fs,mode,cmap,FN,W,H)
[nf,nc] = size(fs);
nn = height(T);
nodes = T;
if ~isempty(fig);fg = figure(fig);else;fg = figure;end
%% convert fs upper triangle weights to matrix form

ix=find(triu(true(nn),1));

if size(fs,2)==1
    adjw = zeros(nn);
    adjg = zeros(nn);
    adjw(ix(fs(:,1)~=0))=fs(fs(:,1)~=0,1);
    adjw = symAdj(adjw,'upper');
    adjg = adjw~=0;
else
    adjg = zeros(nn,nn,nc);
    adjw = zeros(nn,nn,nc);
    [ixI,ixJ]=find(triu(true(nn),1));
    
    for ii=1:nc
        ix = find(fs(:,ii)~=0);
        for jj=1:nnz(ix)
            adjw(ixI(ix(jj)),ixJ(ix(jj)),ii)=fs(ix(jj),ii);
            adjg(ixI(ix(jj)),ixJ(ix(jj)),ii)=1;
        end
    end
end
%% extract filled subnetwork
idx = [];
% pretty disgusting code a head needs revision when i have time
sadjw = adjw;
degree = zeros(nn,1);
weight = zeros(nn,1);
for ii=1:nc
    tmp = symAdj(adjw(:,:,ii),'upper');
    switch mode
        case 'positive';tmp(tmp<0)=0;
        case 'negative';tmp(tmp>0)=0;
    end
    degree = degree+any(tmp,2);% count multiplex degree
    weight = weight+sum(tmp,2);
    [I,J] = find(tmp);
    idx = [idx;[I,J,ones(size(I))*ii]];
end
lidx=sub2ind(size(adjw),idx(:,1),idx(:,2),idx(:,3));

ix = degree>0;
%imagesc(sum(sub_adjw,3))
sub_adj = sadjw(ix,ix,:);% 
%% define order based on AAL code
sT = table(T(ix));
% compute degree based on mode 


nw = sum(sub_adj)';
if nnz(nw~=0)
    clim = round(max(abs(nw)),1).*[-1,1];
    nc = gen.ind2Cmap(cmap,[clim(:);nw]);
    sT.nmap = nc(3:end,:);
end
labels = matlab.lang.makeUniqueStrings(sT.Var1);
%% focus on positive side
%sub_adj(sub_adj<0)= 0;
%num = arrayfun(@(a) sprintf('%02d',a),sT.Var1,'un',0);
if ~exist('cmap','var')||isempty(cmap)
    cmap = hsv(360);
    cmap = genCmap([cmap(220,:);cmap(220,:).^.5;.85,.85,.85;cmap(15,:).^.5;cmap(15,:)]);
end
usub_adj = sub_adjw;
usub_adj(tril(true(size(sub_adj)),-1)) = 0;

A = sparse(usub_adj);
N = size(sub_adj,1);

%% 
clf(fg)
ax = axes('position',[0.25 0.25 0.5 0.5]);
n=size(sub_adj,1);
rad = linspace(0.5*pi,2.5*pi,n+1);
X=cos(rad);
Y=sin(rad);
hold on


X = X(:,1:n);
Y = Y(:,1:n);
t = linspace(0,1,101);

[I,J] = find(A~=0);
w =  c3nl_scale(abs(nonzeros(A)),2,6);
a =  c3nl_scale(abs(nonzeros(A)),0.4,0.65);
clim = round(max(abs(sub_adj(:)))).*[-1,1];
c =  ind2Cmap(cmap,[clim(:);nonzeros(A)]);
c = c(3:end,:);
%% 
for ii=1:nnz(usub_adj)
[x,y] = plot.arc(rad(I(ii)),rad(J(ii)));
 patch([x(:);NaN],[y(:);NaN],'k','linewidth',w(ii),'edgealpha',a(ii),'edgecolor',c(ii,:),'parent', ax);
end

[cx,cy]=circle(0,0,1);% get a unit circle
%plot(cx,cy)
patch([cx';NaN],[cy';NaN],'k','linewidth',3,'edgealpha',0.5,'edgecolor','k','parent', ax);
nw = sum(sub_adj)';
d = sum(sub_adj~=0,2)';
nws =  c3nl_scale(nw,0.02,0.06);
clim = round(max(abs(sub_adj(:))),1).*[-1,1];
nc = ind2Cmap(cmap,[clim(:);nw]);
nc = nc(3:end,:);
for ii=1:n
    xx = X(ii);yy = Y(ii);
    [cx,cy]=plot.circle(xx,yy,nws(ii));% get small circles
    plot(cx,cy,'k','linewidth',2,'parent', ax)
    patch(cx,cy,'k','FaceAlpha',1,'facecolor',nc(ii,:),'parent', ax);% plot the nodes
    t = atan2(yy,xx);
    if abs(t) > pi/2        
        text(1.1*xx,1.1*yy,labels{ii},'rotation', 180*(t/pi + 1),'HorizontalAlignment','right','interpreter','none','parent', ax)
    else
        text(1.1*xx,1.1*yy,labels{ii},'rotation', t*180/pi,'HorizontalAlignment','left','interpreter','none','parent', ax)
    end
end

axis image
axis 'off'

savePDF(FN,W,H);
end