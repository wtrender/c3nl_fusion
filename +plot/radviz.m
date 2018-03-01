function radviz(yc,yp,labels,cmap,gmap,ms,ax)
label = 1;
disp('');

if ~exist('ax','var')
    figure;
    ax = axes('position',[0.25 0.25 0.5 0.5]);
end

labels = labels(:);
if size(yc,2)==1
    g = grp2idx(yc);
    Y = yc;
else 
    g = grp2idx(yc(:,1)); % ground truth
    Y = yc(:,1);
    gp = grp2idx(yc(:,2)); % predicted class
end
nc = numel(categories(Y));
[~,n] = size(yp);
ang=0.5*pi:2*pi/n:2.5*pi;
X=cos(ang);
Y=sin(ang);
[cx,cy]=plot.circle(0,0,1);% get a unit circle
plot(cx,cy,'k-','linewidth',2,'parent', ax);hold on
for ii=1:n
    xx = X(ii);yy = Y(ii);
    [cx,cy]=plot.circle(xx,yy,0.07);% get small circles
    plot(cx,cy,'k','linewidth',2,'parent', ax)
    patch(cx,cy,'k','FaceAlpha',1,'facecolor',cmap(ii,:),'parent', ax);% plot the nodes
    if label
    t = atan2(yy,xx);
    if abs(t) > pi/2        
        text(1.1*xx,1.1*yy,labels{ii},'rotation', 180*(t/pi + 1),'HorizontalAlignment','right','interpreter','none','parent', ax)
    else
        text(1.1*xx,1.1*yy,labels{ii},'rotation', t*180/pi,'HorizontalAlignment','left','interpreter','none','parent', ax)
    end
    end
end



axis([-1.2,1.2,-1.2,1.2]);

x = cos(ang(1:end-1))*yp';
y = sin(ang(1:end-1))*yp';

% y = fft(M);                               % Compute DFT of x
% m = abs(y);                               % Magnitude
% p = unwrap(angle(y));                     % Phase
hold on


if ~exist('gp','var')
    for ii=1:nc
        ix = g==ii;
        plot(x(ix) ,y(ix),'o','markerfacecolor',gmap(ii,:),'markersize',ms,'markeredgecolor','none');
    end
else
   ix = (gp == g);
   for ii=1:nc
        ix1 = ((g==ii).*ix)>0;
        if any(ix1)
            plot(x(ix1) ,y(ix1),'o','markerfacecolor',gmap(ii,:),'markersize',ms,'markeredgecolor','none');
        end
   end
   for ii=1:nc
        ix1 = ((g==ii).*~ix)>0;
        if any(ix1)
            p = plot(x(ix1) ,y(ix1),'*','markerfacecolor',gmap(ii,:),'markersize',ms,'markeredgecolor',gmap(ii,:));
        end
    end 
end

axis off
end
