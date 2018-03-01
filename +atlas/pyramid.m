function P=pyramid(Y,k)
dim = size(Y);
nd = numel(dim);
[~,~,ker]= c3nl.neighb(dim,1);
P = cell(k+1,1);
P{1} = Y;
c = repmat({':'},nd-1,1);
for ii=1:k
   tmp = convn(P{ii},ker,'same');
   for jj=1:nd % reduce by half
       tmp = shiftdim(tmp,1);
       tmp = tmp(1:2:end,c{:});
   end
   P{ii+1} = tmp;
end

end

