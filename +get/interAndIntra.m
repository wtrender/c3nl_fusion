function [AB,tmap] = interAndIntra(n,gp,cmap)

ix1 = find(triu(true(n),1));
[I,J]=find(triu(true(n),1));
LUT = zeros(2^8,1,'uint8');
LUT(2:numel(gp)+1) = gp;
A = zeros(n);
B = zeros(n);
A(ix1) = intlut(uint8(I),LUT);
B(ix1) = intlut(uint8(J),LUT);
A = apply.symAdj(A,'upper');
B = apply.symAdj(B,'upper');
pr = [[1:max(gp);1:max(gp)]';nchoosek(1:max(gp),2)];
AB = zeros(n);
c = 1;
tmap = zeros(size(pr,1),3);
for ii=1:size(pr,1)
    tmp = A==pr(ii,1)&B==pr(ii,2)|A==pr(ii,2)&B==pr(ii,1);
    AB = AB+tmp*c;
    c=c+1;
    if diff(pr(ii,:))
        tmp = get.cmap([cmap(pr(ii,1),:);cmap(pr(ii,2),:)],3);
        tmap(ii,:) = tmp(2,:);
    else
        tmap(ii,:) = cmap(ii,:);
    end
end