function [MU,STD,LEVELS,tmap,AB,ix2] = meanconn(X,gp,n,T,cmap,out)
ix1 = find(triu(true(n),1));
[I,J]=find(triu(true(n),1));
LUT = zeros(2^8,1,'uint8');
LUT(T.number+1) = T.cid;
A = zeros(n);
B = zeros(n);
A(ix1) = intlut(uint8(I),LUT);
B(ix1) = intlut(uint8(J),LUT);
A = apply.symAdj(A,'upper');
B = apply.symAdj(B,'upper');
[~,ix2] = sort(T.cid);
pr = [[1:max(T.cid);1:max(T.cid)]';nchoosek(1:max(T.cid),2)];
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
M = repmat(zeros(n),1,1,size(X,1));
for ii=1:size(X,1)
    tmp = zeros(n);
    tmp(ix1) = X(ii,:);
    M(:,:,ii)= apply.symAdj(tmp,'upper');
end
Mr = reshape(M,[],size(X,1));
ABr = reshape(AB,[],1);
GP = dummyvar(gp);
MU = zeros(max(ABr),size(GP,2));
STD = zeros(max(ABr),size(GP,2));
for ii=1:max(ABr)
    ix = ABr==ii;
    MU(ii,:)=mean(Mr(ix,:)*GP./sum(GP));
    STD(ii,:)=std(Mr(ix,:)*GP./sum(GP));
end
LEVELS = categories(gp);

end

