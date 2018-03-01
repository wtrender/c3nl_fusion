function L = gbs(Y,k)
%% ATLAS.GBS: One line description of what the function or script performs
%
%   __           _             
%  / _|         (_)            
% | |_ _   _ ___ _  ___  _ __    
% |  _| | | / __| |/ _ \| `_ \    :- Functional and Structural 
% | | | |_| \__ \ | (_) | | | |      Integration of Neuroimages
% |_|  \__,_|___/_|\___/|_| |_|
%
%
%% AUTHOR:  Eyal Soreq
%  EMAIL:  e.soreq14@imperial.ac.uk
%  AFFILIATION:  Imperial College London
%  VERSION:  0.0 CREATED: 10-Oct-2017 08:31:57
%
%% INPUTS:
%    input01 - 
%    input02 - 
%    input03 - 
%    input04 - 
%
%
%% OUTPUT:
%
%% EXAMPLES:
%{
gbs(input01,input02,input03,input04)
%}
%
%% DEPENDENCIES:
%
% This file is part of Fusion Pipeline
% Fusion Pipeline is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% Fusion Pipeline is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with Fusion Pipeline.If not, see <http://www.gnu.org/licenses/>.
%------------- BEGIN CODE --------------
%
disp('')
Y = double(Y);
o_dim = size(Y);
nd = numel(o_dim);
Y = padarray(Y,ones(numel(o_dim),1),nan); % pad image to account for bounderies
Pdim =size(Y);
Y(Y<0.4)=nan;
d = c3nl.neighb(Pdim,1,'city');
ix = find(~isnan(Y));
nv = numel(ix);
nc = nv;
L = zeros(Pdim);
t = zeros(Pdim);
t(ix) = k/1;

L(ix) = 1:nv;
ew = zeros(numel(ix),numel(d));
for ii=1:nv
  ew(ii,:)= abs(Y(ix(ii))-Y(ix(ii)+d)); 
end
[~,idx] = sort(ew(:));
tic
for ii=1:numel(idx)
   if ~isnan(ew(idx(ii)))
      [ixr,ixc]=ind2sub(size(ew),idx(ii));
      a = ix(ixr);
      b = ix(ixr)+d(ixc);
      if (L(a)~=L(b))
         if ew(idx(ii))<=t(a) && ew(idx(ii))<=t(b)
            L(L==L(b))=L(a);
            t(L==L(b))=ew(idx(ii))+k/nnz(L==L(b));
         end
      end
    if ~mod(ii,100);fprintf('*');end
    if ~mod(ii,6000);fprintf('\n');end
   end

end
toc


c = repmat({':'},nd-1,1);
for ii=1:nd % crop image back to normal
    L = shiftdim(L,1);
    L = L(2:end-1,c{:});
end
% figure(1);plot.imOverlay(get.montage(Y,3),get.montage(L,3),0.5)

end
%------------- END OF CODE --------------
