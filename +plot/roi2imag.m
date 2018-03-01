function Output = roi2imag(D,obj)

[type,pw,ph,az,el,cmap,G,Alpha,Smooth,cortex,lim,light,T,adj,r,thr,marker,bkg] = getParams(obj);

fg = figure(99);clf;
set(fg,'units','pixels','position',[0 0 pw ph],'Visible','on','resize','off')
%fg = figure('Number',1,'units','pixels','position',[0 0 pw ph],'Visible','on','resize','off');
ax=axes('position' ,[0.05 0.05 .9 .9]);axis off;
switch type
     case {'adj'}
        % cortex = the brain shell
        % net  = the adjency metrix
        % T = the table of info relating to the net
        %adj = symAdj(adj,'upper');
        if ~nnz(adj(:));Output=[];return;end
        ix = zeros(height(T),1);cent=[];
        %% plot tubes
        patch('Vertices',cortex.V,'Faces',cortex.F,'FaceColor','r','FaceAlpha',0,'EdgeColor','none');
        %        figure;ax=axes('position' ,[0.05 0.05 0.9 0.9]);axis off;
        hold on
        xyz =     T.Centroid;
        sz = size(adj);
        ix1 = adj ~= 0 &~isnan(adj); w = zeros(sz);
        Lmap = zeros(numel(w),3); % generate link color map
        if ~exist('cmap','var')
            cmap = hsv(360);
            o2b = genCmap([cmap(220,:);cmap(220,:).^.5;.85,.85,.85;cmap(15,:).^.5;cmap(15,:)]);
        else 
            o2b = cmap;
        end
        %o2b = flipud([244 122 32;247 146 71;250 171 113;252 197 154;250 250 250;148 187 227;98 161 215;8 140 205;0 123 195]);
        %p2g = [237 30 141;237 78 156;236 116 173;244 155 194;250 250 250;162 210 146;124 195 98;83 184 72;66 182 73];
        clim = round(max(abs(adj(:))),1).*[-1,1];
        gmap =[];
        if ~strcmpi(type,'adj_alpha')
            tmp =  ind2Cmap(o2b,[clim(:);adj(ix1)]);
            Lmap(ix1,:) = tmp(3:end,:);
            Lmap = reshape(Lmap,[sz,3]);
        else
            Lmap = ones([sz,3]);
        end
        w(ix1) = c3nl_scale(abs(adj(ix1)),.5,2.5);
        for ii=1:size(adj,2) % go over the nodes
            if any(abs(adj(ii,:))>0)
                ix1 = find(abs(adj(ii,:))>0);
                for jj=1:length(ix1)
                    tmp = xyz([ii,ix1(jj)],[2,1,3])';
                    c = squeeze(Lmap(ii,ix1(jj),:));
                    [X,Y,Z] = genTube(tmp,15,w(ii,ix1(jj)));
                    surf(ax,X,Y,Z,'FaceLighting','gouraud','facecolor',c,'edgecolor','none')
                end
            end
        end
        tmp = sum(sign(abs(symAdj(adj,'upper'))),2);
        tmp(tmp~=0) = c3nl_scale(tmp(tmp~=0),2.5,6.5);
        for ii=1:size(adj,2) % go over the nodes
            
            if tmp(ii)>0 % if it has any connections plot the node
                if strcmpi(type,'adj_alpha');cmap = ones(1,3); else cmap = T.nmap(ii,:); end
                [X,Y,Z] = ellipsoid(xyz(ii,2),xyz(ii,1),xyz(ii,3),tmp(ii),tmp(ii),tmp(ii),20);
                surf(ax,X,Y,Z,'FaceLighting','gouraud','facecolor',cmap,'faceAlpha',1,'edgecolor','none')
                cent = [cent;xyz(ii,:)];
                ix(ii) = 1;
                gmap = [gmap;cmap];
                %plotSpheres(xyz(ii,:),2,20,[1,0,0],ax)
%             else 
%                 
%                 disp 'void'
            
            end
            
        end
        view(az,el);axis 'off';daspect([1 1 1]);axis(lim)
        camlight(-80,-10); camlight(80,-10);
        ax.YDir = 'reverse';
        %camproj('perspective')
        tmp = getframe(ax);
        Output.img = tmp.cdata;
        Output.lim = [xlim,ylim,zlim];
        Output.alpha = rgb2gray(tmp.cdata)~=240;
    
    case 'stack'
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
        mmap = cell(size(D,4),1);
        
        switch or
            case 1
                for jj=1:size(D,4)
                    vmx = squeeze(max(D(:,:,:,jj),[],or));
                    tmp = t{or};
                    if size(D,4)>1
                        mc = max(0,cmap(jj,:)-[0.1,0.05,0.5]);
                        tmap = [mc.^10;cmap(jj,:);(mc+0.1).^0.01];                   
                    else 
                        tmap = cmap;
                    end
                    mmap{jj} = tmap;
                    for ii=1:numel(t{or})
                        v = squeeze(D(ii,:,:,jj));
                        if numel(v)~=sum(isnan(v(:)))
                            v= vmx.*(v>0);
                            vc = [v(~isnan(v));max(D(:));min(D(:))];
                            if numel(unique(vc))>2
                                c = nan([numel(v),3]);
                                cc = ind2Cmap(tmap,vc);
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
                            vc = [v(~isnan(v));max(D(:));min(D(:))];
                            if numel(unique(vc))>2
                                c = nan([numel(v),3]);
                                cc = ind2Cmap(tmap,vc);
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
                            vc = [v(~isnan(v));max(D(:));min(D(:))];
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
        tmp = getframe(ax);
        Output.img = tmp.cdata;%imgaussfilt(tmp.cdata,2);
        Output.lim = [xlim,ylim,zlim];
        Output.alpha = rgb2gray(tmp.cdata)~=240;
        Output.mmap = mmap;
    case 'mip'
        %
        %         case 'axial';az = 90*dir;el=90;
        %         case 'sagittal';if (dir==3); az = 180;else; az = 360;end; el=360;
        %         case 'coronal';if (dir==3); az = 90;else; az = 270;end; el=360;
        y=round(sin(az*pi/180)*cos(el*pi/180),2);
        x=round(cos(az*pi/180)*cos(el*pi/180),2);
        z=round(sin(el*pi/180),2);
        dim = find([x,y,z]);
        D(D<=thr) = NaN;
        v = squeeze(max(D,[],dim));
        % D(isnan(D)) = 0;
        % v = squeeze(sum(D,dim));
        [~,b]=ismember(size(v),size(D));
        t{1} = G.xmm;t{2} = G.ymm;t{3} = G.zmm;
        [X,Y] = meshgrid(t{b(2)},t{b(1)});
        switch dim
            case 1; surf(Y,zeros(size(X)),X,v,'edgecolor','none','parent',ax);
            case 2; surf(zeros(size(X)),Y,X,v,'edgecolor','none','parent',ax);ax.YDir = 'reverse';
            case 3; surf(X,Y,zeros(size(X)),v,'edgecolor','none','parent',ax);ax.YDir = 'reverse';
        end
        
        view(az,el);axis 'off';daspect([1 1 1]);axis(lim)
        fg.Color = [0,0,0];
        colormap(cmap)
        tmp = getframe(ax);
        Output.img = tmp.cdata;%imgaussfilt(tmp.cdata,2);
        Output.lim = [xlim,ylim,zlim];
        Output.alpha = rgb2gray(tmp.cdata)>0;
        
    case 'pc'
        mypcshow(ax,D,G,thr,marker,cmap,bkg);
        view(az,el);axis 'off';daspect([1 1 1]);
        if ~isempty(lim);axis(lim);end
        tmp = getframe(ax);
        Output.alpha = rgb2gray(tmp.cdata)>0;%imgaussfilt((rgb2gray(tmp.cdata)>0)*1.0,2);
        if ~isempty(lim);axis(lim);end
        Output.img = tmp.cdata;%imgaussfilt(tmp.cdata,2);
        Output.lim = [xlim,ylim,zlim];
        
    case {'AO','AO_alpha'}
        co=[1,1,1];ss=.5;ds=.1;as=.8;
        if strcmpi(type,'AO'); fvc = bsxfun(@times,repmat(co,size(cortex.V,1),1), 1-cortex.AO);
        else fvc = bsxfun(@times,repmat(co,size(cortex.V,1),1), cortex.AO); end
        patch('Vertices',cortex.V,'Faces',cortex.F,'FaceVertexCData',fvc,'FaceColor','interp','EdgeColor','none','SpecularStrength',ss,'DiffuseStrength',ds,'AmbientStrength',as,'parent',ax);
        view(az,el);axis 'off';daspect([1 1 1]);
        if ~isempty(lim);axis(lim);end
        if light;camlight(-80,-10); camlight(80,-10);camlight;lighting phong;end
        tmp = getframe(ax);
        Output.img = c3nl_scale(double(rgb2gray(tmp.cdata)));
        Output.lim = [xlim,ylim,zlim];
    case 'Alpha'
        patch('Vertices',cortex.V,'Faces',cortex.F,'FaceColor','w','EdgeColor','none');
        view(az,el);axis 'off';daspect([1 1 1]);
        if ~isempty(lim);axis(lim);end
        tmp = getframe(ax);
        Output = c3nl_scale(double(rgb2gray(tmp.cdata)));
    case 'spec'
        patch('Vertices',cortex.V,'Faces',cortex.F,'FaceColor','k','EdgeColor','none','facealpha',1,'facelighting','phong','SpecularColorReflectance',1 ,'DiffuseStrength',0.3,'SpecularExponent',10,'SpecularStrength',1);
        view(az,el);axis 'off';daspect([1 1 1]);
        if light;camlight(-80,-10); camlight(80,-10);camlight;lighting phong;end
        if ~isempty(lim);axis(lim);end
        tmp = getframe(ax);
        Output.img = rgb2gray(tmp.cdata);
        Output.lim = [xlim,ylim,zlim];
    case {'roi','roi_alpha'}
        if strcmpi(type,'Alpha');cmap = ones(length(unique(D)),3);end
        collateSurface(D,G,Alpha,Smooth,cmap);
        view(az,el);axis 'off';
        if ~isempty(lim);axis(lim);end
        daspect([1 1 1]);ax.YDir = 'reverse';
        if strcmpi(type,'roi_alpha');tmp = getframe(ax);Output.img = c3nl_scale(double(rgb2gray(tmp.cdata)));Output.lim = [xlim,ylim,zlim];
        else
            camlight;lighting phong;
            tmp = getframe(ax);
            Output.img = tmp.cdata;
            Output.lim = [xlim,ylim,zlim];
            Output.alpha = rgb2gray(tmp.cdata)~=240;
            
        end
    case 'flat_roi'
        collateSurface(D,G,Alpha,Smooth,cmap);
        view(az,el);axis 'off';
        if ~isempty(lim);axis(lim);end
        daspect([1 1 1]);ax.YDir = 'reverse';
        tmp = getframe(ax);
        Output.img = tmp.cdata;
        Output.lim = [xlim,ylim,zlim];
        Output.alpha = rgb2gray(tmp.cdata)~=240;
            
    case {'nodes','nodes_alpha'}
        xyz =     T.Centroid;
        hold on
        for ii=1:size(adj,2) % go over the nodes
            if any(abs(adj(ii,:))>0) % if it has any connections plot the node
                if strcmpi(type,'nodes_alpha');map = ones(1,3); else map = cmap(ii,:); end
                [X,Y,Z] = ellipsoid(xyz(ii,2),xyz(ii,1),xyz(ii,3),r,r,r,20);
                surf(ax,X,Y,Z,'FaceLighting','gouraud','facecolor',map,'faceAlpha',1,'edgecolor','none')
            end
        end
        hold off
        view(az,el);axis 'off';ax.YDir = 'reverse';
        if ~isempty(lim);axis(lim);end
        daspect([1 1 1]);
        if strcmpi(type,'nodes_alpha');tmp = getframe(ax);Output.img = c3nl_scale(double(rgb2gray(tmp.cdata)));Output.lim = [xlim,ylim,zlim];
        else
            camlight(-80,-10); camlight(80,-10);camlight;lighting phong;
            tmp = getframe(ax);
            Output.img = tmp.cdata;
            Output.lim = [xlim,ylim,zlim];
            Output.alpha = rgb2gray(tmp.cdata)~=240;
        end
        
end
close(fg)
end



function [type,pw,ph,az,el,cmap,G,Alpha,Smooth,cortex,lim,light,T,adj,r,thr,marker,bkg] = getParams(obj)

prop = {'light','lim','type','dpi','cw','ch','w','h','pw','ph','az','el','cmap','cortex','G','direction','view','Alpha','Smooth','T','adj','r','thr','marker','bkg'};
if ~isstruct(obj);ia = zeros(numel(prop),1); else [ia,~] = ismember(prop,fieldnames(obj));end
for ii=1:length(ia)
    tmp = prop{ii};
    switch tmp
        case 'r';if ia(ii);r=obj.(tmp);else; r=4;end % res in dpi default = 96
        case 'T';if ia(ii);T=obj.(tmp);else; T=[];end % res in dpi default = 96
        case 'adj';if ia(ii);adj=obj.(tmp);else; adj=[];end % res in dpi default = 96
        case 'light';if ia(ii);light=obj.(tmp);else; light=1;end % res in dpi default = 96
        case 'lim';if ia(ii);lim=obj.(tmp);else; lim=[-120,80,-80,80,-80,100];end % res in dpi default = 96
        case 'Smooth';if ia(ii);Smooth=obj.(tmp);else; Smooth=1;end % res in dpi default = 96
        case 'Alpha';if ia(ii);Alpha=obj.(tmp);else; Alpha=0.5;end % res in dpi default = 96
        case 'cortex';if ia(ii);cortex=obj.(tmp);else; load ('+atlas/cortex');end % res in dpi default = 96
        case 'type';if ia(ii);type=obj.(tmp);else; type='roi';end % res in dpi default = 96
        case 'cmap';if ia(ii);cmap=obj.(tmp);else; cmap=jet;end % res in dpi default = 96
        case 'G';if ia(ii);G=obj.(tmp);else; G=cortex.G;end % todo create a deafult grid settings
        case 'dpi';if ia(ii);dpi=obj.(tmp);else; dpi=96;end % res in dpi default = 96
        case 'cw';if ia(ii);cw=obj.(tmp);else; cw=20;end % defult res is 20 cm
        case 'ch';if ia(ii);ch=obj.(tmp);else; ch=20;end % defult res is 20 cm
        case 'w';if ia(ii);w=obj.(tmp);else; w=1;end %
        case 'h';if ia(ii);h=obj.(tmp);else; h=w;end %
        case 'ph';if ia(ii);ph=obj.(tmp);else; ph=fix(dpi*ch*h/2.54);end %
        case 'pw';if ia(ii);pw=obj.(tmp);else; pw=fix(dpi*cw*w/2.54);end %
        case 'az';if ia(ii);az=obj.(tmp);else; az=270;end
        case 'el';if ia(ii);el=obj.(tmp);else; el=90;end
        case 'direction';if ia(ii);dir=obj.(tmp);else dir=3;end
        case 'thr';if ia(ii);thr=obj.(tmp);else; thr=3;end
        case 'marker';if ia(ii);marker=obj.(tmp);else marker='*';end
        case 'bkg';if ia(ii);bkg=obj.(tmp);else; bkg=[0,0,0];end
            
        case 'view'
            if ia(ii)
                switch obj.(tmp)
                    case 'axial';az = 90*dir;el=90;
                    case 'sagittal';if (dir==3); az = 180;else; az = 360;end; el=360;
                    case 'coronal';if (dir==3); az = 90;else; az = 270;end; el=360;
                end
            end
    end
end
end