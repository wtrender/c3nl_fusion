function y = net(d,Y,dist)
    y = Y(1:dist:d(1),1:dist:d(2),1:dist:d(3));
%     L = atlas.watershed(y,1,5);
%     plot.montageOverlay(y,L,0.5)
%     numel(unique(L))
%     
%     stats = struct2array(regionprops(L,{'area'}))';
%     ix = find(stats<4);
%     plot.montageOverlay(y,ismember(L,ix),0.3)

            
end