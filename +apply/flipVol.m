function [V,G]= flipVol(v,g,tflip)


G = g;
V = v;
if ~isempty(tflip)
    for ii=1:numel(tflip)
        switch tflip(ii)
            case 1 
                G.mm{1} = fliplr(G.mm{1});
                V = flip(V,1);
                [~,ix]=max(abs(G.mat(1,1:3)));                
                G.vx(1) = -G.vx(1);
                G.mat(1,[ix,4]) =   [-1*G.mat(ix,1),G.mm{1}(1)-G.vx(1)];
                G.bb(1,1:2)=G.bb(1,[2,1]);
            case 2
                G.mm{2} = fliplr(G.mm{2});
                V = flip(V,2);
                [~,ix]=max(abs(G.mat(2,1:3)));                
                G.vx(2) = -G.vx(2);
                G.mat(2,[ix,4]) =   [-1*G.mat(ix,2),G.mm{2}(1)-G.vx(2)];
                G.bb(2,1:2)=G.bb(2,[2,1]);
            case 3
                G.mm{3} = fliplr(G.mm{3});
                V = flip(V,3);
                [~,ix]=max(abs(G.mat(3,1:3)));                
                G.vx(3) = -G.vx(3);
                G.mat(3,[ix,4]) =   [-1*G.mat(ix,3),G.mm{3}(1)-G.vx(3)];
                G.bb(3,1:2)=G.bb(3,[2,1]);
        end
    end
end



end