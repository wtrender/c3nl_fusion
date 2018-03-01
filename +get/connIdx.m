function [AB,ix2,pr] = connIdx(n,cid)
%% GET.CONNIDX: One line description of what the function or script performs
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
%  VERSION:  0.0 CREATED: 29-Nov-2017 17:19:44
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
connIdx(input01,input02,input03,input04)
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
ix1 = find(triu(true(n),1));
[I,J]=find(triu(true(n),1));
LUT = zeros(2^8,1,'uint8');
LUT((1:numel(cid))+1) = cid;
A = zeros(n);
B = zeros(n);
A(ix1) = intlut(uint8(I),LUT);
B(ix1) = intlut(uint8(J),LUT);
A = apply.symAdj(A,'upper');
B = apply.symAdj(B,'upper');
[~,ix2] = sort(cid);
pr = [[1:max(cid);1:max(cid)]';nchoosek(1:max(cid),2)];
AB = zeros(n);
pr(:,3) = (1:size(pr,1))';
for ii=1:size(pr,1)
    tmp = A==pr(ii,1)&B==pr(ii,2)|A==pr(ii,2)&B==pr(ii,1);
    AB = AB+tmp*pr(ii,3);
end
end
%------------- END OF CODE --------------
