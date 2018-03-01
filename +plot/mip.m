function [Output,mmap] = mip(D,G,dim,thr,cmap)
switch dim
    case 'a';az = 270;el=90;%axial
    case 'rs';az = 360; el=360;%sagittal
    case 'ls';az = 180; el=360;%sagittal
    case 'c';az = 270; el=360;%coronal
end
        fg = figure(99);clf;
        set(fg,'units','pixels','position',[0 0 600 600],'Visible','on','resize','off');
        ax=axes('position' ,[0 0 1 1]);axis off;
        load +atlas/cortex
        patch('Vertices',cortex.V,'Faces',cortex.F,'FaceVertexCData',1-cortex.AO,'FaceColor','interp','EdgeColor','none','FaceAlpha',1);
        view(az,el);axis 'off';daspect([1 1 1]);colormap gray
        light
        tmp = getframe(fg);
        Output.cortex = tmp.cdata;
        lim = axis;
        cla(ax)
        y=round(sin(az*pi/180)*cos(el*pi/180),2);
        x=round(cos(az*pi/180)*cos(el*pi/180),2);
        z=round(sin(el*pi/180),2);
        dim =[x,y,z];
        or = find(dim);
        D(D<=thr) = NaN;
        b = find(~dim);
        t{1} = G.mm{1};t{2} = G.mm{2};t{3} = G.mm{3};
        [X,Y] = meshgrid(t{b(2)},t{b(1)});
        
        hold(ax,'on')
        mx = max(D(:));
        mn = min(D(:));
        mmap = cell(size(D,4),1);
        
        switch or
            case 1
                for jj=1:size(D,4)
                    vmx = squeeze(max(D(:,:,:,jj),[],or));
                    tmp = t{or};
                    if size(D,4)>1
                        mc = min(max(jet(12)-0.1,0.01),0.95);
                        tmap = c3nl.cmap([mc.^20;mc;mc.^.01],100);
                    else 
                        tmap = cmap;
                    end
                    mmap{jj} = tmap;
                    for ii=1:numel(t{or})
                        v = squeeze(D(ii,:,:,jj));
                        if numel(v)~=sum(isnan(v(:)))
                            v= vmx.*(v>0);
                            vc = [v(~isnan(v));mx+0.01*mn;mn];
                            if numel(unique(vc))>2
                                c = nan([numel(v),3]);
                                cc = c3nl.ind2Cmap(tmap,vc);
                                ix = find(~isnan(v));
                                c(ix,:) = cc(1:end-2,:);
                                v(isnan(v))=0;
                                surf(Y,ones(size(X))*tmp(ii)+(0.1*ii),X,v,'AlphaData',v.^.2,'Cdata',reshape(c,[size(v),3]),'FaceAlpha','flat','edgecolor','none','parent',ax);
                                %surf(Y,ones(size(X))*tmp(ii)+(0.1*ii),X,v,'Cdata',reshape(c,[size(v),3]),'edgecolor','none','parent',ax);
                            end
                        end
                    end
                end
                
                
            case 2
                for jj=1:size(D,4)
                    tmp = t{or};
                    vmx = squeeze(max(D(:,:,:,jj),[],or));
                    if size(D,4)>1
                        mc = max(0,cmap(jj,:)-[0.1,0.05,0.5]);
                        tmap = [mc.^10;cmap(jj,:);(mc+0.1).^0.01];                   
                    else 
                        tmap = cmap;
                    end
                    mmap{jj} = tmap;
                    for ii=1:numel(t{or})
                        v = squeeze(D(:,ii,:,jj));
                        if numel(v)~=sum(isnan(v(:)))
                            v= vmx.*(v>0);
                            vc = [v(~isnan(v));;mx;mn];
                            if numel(unique(vc))>2
                                c = nan([numel(v),3]);
                                cc = c3nl.ind2Cmap(tmap,vc);
                                ix = find(~isnan(v));
                                c(ix,:) = cc(1:end-2,:);
                                v(isnan(v))=0;
                                surf(ones(size(X))*tmp(ii),Y,X,v,'AlphaData',v.^.2,'Cdata',reshape(c,[size(v),3]),'FaceAlpha','flat','edgecolor','none','parent',ax);ax.YDir = 'reverse';%flat
                                %surf(ones(size(X))*tmp(ii),Y,X,v,'Cdata',reshape(c,[size(v),3]),'edgecolor','none','parent',ax);ax.YDir = 'reverse';%flat
                            end
                        end
                    end
                end
            case 3
                for jj=1:size(D,4)
                    tmp = t{or};
                    vmx = squeeze(max(D(:,:,:,jj),[],or));
                     if size(D,4)>1
                        mc = max(0,cmap(jj,:)-[0.1,0.05,0.5]);
                        tmap = [mc.^10;cmap(jj,:);(mc+0.1).^0.01];                   
                    else 
                        tmap = cmap;
                    end
                    mmap{jj} = tmap;
                    for ii=1:numel(t{or})
                        v = squeeze(D(:,:,ii,jj));
                        if numel(v)~=sum(isnan(v(:)))
                            v= vmx.*(v>0);
                            vc = [v(~isnan(v));mx;mn];
                            if numel(unique(vc))>2
                                c = nan([numel(v),3]);
                                cc = c3nl.ind2Cmap(tmap,vc);
                                ix = find(~isnan(v));
                                c(ix,:) = cc(1:end-2,:);
                                v(isnan(v))=0;

                                surf(X,Y,ones(size(X))*tmp(ii),v,'AlphaData',v.^.2,'Cdata',reshape(c,[size(v),3]),'FaceAlpha','flat','edgecolor','none','parent',ax);ax.YDir = 'reverse';
                                %surf(X,Y,ones(size(X))*tmp(ii),v,'Cdata',reshape(c,[size(v),3]),'edgecolor','none','parent',ax);ax.YDir = 'reverse';
                                
                            end
                        end
                    end
                end
        end
        view(az,el);axis 'off';daspect([1 1 1]);axis(lim)
        %fg.Color = [0,0,0];
        tmp = getframe(fg);
        Output.img = tmp.cdata;%imgaussfilt(tmp.cdata,2);
        Output.lim = [xlim,ylim,zlim];
        Output.alpha = rgb2gray(tmp.cdata)~=240;
        Output.mmap = mmap;
        cla(ax)
        patch('Vertices',cortex.V,'Faces',cortex.F,'FaceColor','w','EdgeColor','none','FaceAlpha',1);
        view(az,el);axis 'off';daspect([1 1 1]);axis(lim);colormap gray
        tmp = getframe(fg);
        Output.mask = (rgb2gray(tmp.cdata)~=240).*1;
        patch('Vertices',cortex.V,'Faces',cortex.F,'FaceVertexCData',cortex.AO,'FaceColor','interp','EdgeColor','none','FaceAlpha',1);
        tmp = getframe(fg);
        Output.Ao = rgb2gray(tmp.cdata);
end