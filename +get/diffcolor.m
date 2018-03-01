function cmap = diffcolor(h,s,v,k,fig)
% cmap = helper.diffcolor([20,240],[0.5,0.8],[0.6,0.9],5,1)
if ~exist('h','var');h=[0,360];end % hue range
if ~exist('s','var');s=[.2,1];end % saturation range we are forcing some color here
if ~exist('s','var');v=[0,1];end % value range
if ~exist('k','var');k=10;end % number of colours

xh = linspace(h(1),h(2),k)/360;
xs = apply.scale(sin(linspace(-(k)*pi,(k)*pi,k)),s(1),s(2));
xv = apply.scale(-cos(linspace(-(k)*pi,(k)*pi,k)),v(1),v(2));
cmap = hsv2rgb([xh;xs;xv]');
if exist('fig','var')
    figure(fig);clf;
    imagesc(1:k);colormap(cmap)
end
% P=sobolset(2);% get quasirandom number generator
% P.Leap = 0;
% 
% X = net(P,k*2);% sample k*2 colors
% %X = rand(k*2,2);
% X(:,1) = c3nl_scale(1:k*2,h(1)/360,h(2)/360); % constrain hue to requested range
% X(:,2) = c3nl_scale(X(:,1),s(1),s(2)); % constrain saturation to requested range
% X(:,3) = c3nl_scale(X(:,2),v(1),v(2)); % constrain value to requested range
% D = pdist(X);
% Z = linkage(D,'weighted');%weighted,average,complete ,'centroid','median','single','ward'
% %[~,id]=max(Y(:,2))
% ord = optimalleaforder(Z,D);
% x = X(ord(1:2:k*2),:);
% cmap = hsv2rgb(x);
% gmap = rgb2gray(cmap);
% x(:,3)=c3nl_scale(x(:,3)./gmap(:,1),v(1),v(2));
% amap=hsv2rgb(x);
% 
% ms = sin(linspace(-(k/2)*pi,(k/2)*pi,k));
% x(:,3)=c3nl_scale(x(:,3).*ms',v(1),v(2));
% 
% mmap=hsv2rgb(x);
% sm = cos(linspace(-pi,pi,k));
% x(:,2)=c3nl_scale(x(:,2).*sm',s(1),s(2));
% msmap=hsv2rgb(x);
% if exist('fig','var')
%     figure(fig)
%     clf
%     hmap = hsv(360);
%     ax1 = subplot(6,1,1);
%     imagesc(1:k);colormap(ax1,hmap(round(linspace(h(1),h(2),k)),:))
%     ax2 = subplot(6,1,2);
%     imagesc(1:k);colormap(ax2,cmap)
%     ax3 = subplot(6,1,3);
%     imagesc(1:k);colormap(ax3,rgb2gray(cmap))
%     ax4 = subplot(6,1,4);
%     imagesc(1:k);colormap(ax4,amap)
%     ax5 = subplot(6,1,5);
%     imagesc(1:k);colormap(ax5,mmap)
%     ax6 = subplot(6,1,6);
%     imagesc(1:k);colormap(ax6,msmap)
% end

