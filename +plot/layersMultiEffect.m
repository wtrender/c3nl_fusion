function img = layersMultiEffect(vol,G,thr,fig,cmap,labels,filename,width, height,cr)
if ~exist('cr','var');cr=0;end
if ~exist('thr','var');thr=2;end
if ~exist('T','var');T=[];end
if ~exist('fig','var');fig=1;end
if ~exist('method','var');method='mip';end
if ~exist('bkg','var');bkg=[246,246,256]./256;end
d = size(vol);

[im.p4,mmap] = plot.layersGlass(vol,G,'a',thr,cmap,600,600,[],[],cr);
tmp = vol;tmp(round(d(1)/2):end,:,:,:) = 0;id = unique(tmp);
im.p1 = plot.layersGlass(tmp,G,'rs',thr,cmap,600,600,[],[],cr);
tmp = vol;tmp(1:round(d(1)/2),:,:,:) = 0;id = unique(tmp);
im.p2 = plot.layersGlass(tmp,G,'ls',thr,cmap,600,600,[],[],cr);
im.p3 = plot.layersGlass(vol,G,'c',thr,cmap,600,600,[],[],cr);
im.c = zeros(d,'uint8');
[x,y,w,h] = get.fancyGrid([2;2;size(vol,4)],0.02,0.03,'weighted',[0.4;0.5;0.1]);

fg = figure(fig);clf
f = -0.05;
fc = 'k';
fs = 0.5;
clf
for ii=1:numel(x)
    tmp = ['p',num2str(ii)];
    if ii<=4
        axf.(tmp) = axes('position',[x(ii),y(ii),w(ii),h(ii)]);
        mask = im.(['p',num2str(ii)]).mask  ;
        ixc = any(mask);
        ixr = any(mask,2);
        imagesc(axf.(tmp),im.(tmp).im(ixr,ixc,:),'AlphaData',mask(ixr,ixc));
        xx = xlim(axf.(tmp));yy=ylim(axf.(tmp));
        xx = [xx(1) mean(xx) xx(2)]+[-f*mean(xx),0,f*mean(xx)];
        yy = [yy(1) mean(yy) yy(2)]+[-f*mean(yy),0,f*mean(yy)];
        xx = [xx xx(2)];yy = [yy(2) yy];
        axis(axf.(tmp), 'off');axis(axf.(tmp),'image');
        hold on;
        for jj=1:4
            [xp,yp] = plot.circle(xx(jj),yy(jj),round(min(range(xx),range(yy))*0.05));
            patch(xp,yp,'w','parent',axf.(tmp),'facealpha',0.5);
        end
        
        switch ii
            case 1;text(axf.(tmp),xx,yy,{'P','S','A','I'},'color',fc,'FontUnits', 'normalized','fontsize',fs/8,'Fontweight','bold','HorizontalAlignment','center');
            case 2;text(axf.(tmp),xx,yy,{'A','S','P','I'},'color',fc,'FontUnits', 'normalized','fontsize',fs/8,'Fontweight','bold','HorizontalAlignment','center');
            case 3;text(axf.(tmp),xx,yy,{'L','S','R','I'},'color',fc,'FontUnits', 'normalized','fontsize',fs/8,'Fontweight','bold','HorizontalAlignment','center');
            case 4;text(axf.(tmp),xx,yy,{'L','A','R','P'},'color',fc,'FontUnits', 'normalized','fontsize',fs/8,'Fontweight','bold','HorizontalAlignment','center');
        end
        hold off
    else
        axf.(tmp) = axes('position',[x(ii),y(ii)*1.02,w(ii)*0.9,h(ii)/2.5]);
        rng = round(linspace(thr,max(reshape(vol(:,:,:,ii-4),[],1)),50),1);
        imagesc(axf.(tmp),rng);
        colormap(axf.(tmp),get.cmap(mmap{ii-4}));
        set(axf.(tmp),'YDir','normal','GridColor',fc,'YColor',fc,'YTick','','XTick',[1,50],'FontUnits', 'normalized','fontsize',fs);
        axf.(tmp).XTickLabel = (cellfun(@(x) num2str(x) ,num2cell(rng(axf.(tmp).XTick)),'un',0));
        xx=xlim; yy=ylim;
        text(axf.(tmp),range(xx)/2,1.55,labels{ii-4},'color',fc,'VerticalAlignment','bottom','HorizontalAlignment','center','FontWeight','bold','FontUnits', 'normalized','fontsize',fs*2);
        axf.(tmp).FontSize = 0.7;
        text(axf.(tmp),3,mean(yy),'T','color','w','FontUnits', 'normalized','fontsize',fs,'FontWeight','bold');
    end
end

if exist('filename','var')
    save.pdf(filename,width, height);
end
if nargout==1
    tmp = getframe(gcf);
    img = tmp.cdata;
end
close (fg);
end