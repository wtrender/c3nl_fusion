function A = scale(A,mn,mx,type)
%% C3NL.SCALE: One line description of what the function or script performs
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
%  VERSION:  0.0 CREATED: 09-Jun-2017 17:18:59
%
%% INPUTS:
%    A -  the matrix to scale 
%    mn - minimum defualt = 0
%    mx - maximum defualt = 1
%    type - row,col or full default='full'
%
%
%% OUTPUT:
%    A = the scaled matrix   
%% EXAMPLES:
%{
As = c3nl.scale(A);
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
sz = size(A);
if ~exist('mn','var');mn=0;end
if ~exist('mx','var');mx=1;end
if mx<mn; throw(MException('C3NL:minLargerThanMax','This function requires that min values will be smaller then max'));end
if ~exist('type','var');type = 'full';end
switch type
    case 'row';A = Scale(A',mn,mx);
    case 'full';A = Scale(A(:),mn,mx);
    case 'col'
        tmp = zeros(sz);
        for ii=1:sz(2)
            tmp(:,ii) = Scale(A(:,ii),mn,mx);
        end
        A = tmp;
end

if strcmpi(type,'full');A = reshape(A,sz);end
if strcmpi(type,'row');A = A';end

end


function d = Scale(D,mn,mx)
if any(isnan(D))
    ix = isnan(D);
    d = D(~ix);
else
    d = D;
end
if sum(abs(diff(d)))
    mi = min(d);
    ma = max(d);
    r = range([mi;ma]);
    d = bsxfun(@rdivide, bsxfun(@minus,d,mi),r);
    d = d*range([mx,mn])+mn;
else 
    d= d*mn;
end
if exist('ix','var')
   tmp = NaN(length(D),1);
   tmp(~ix) = d;
   d= tmp;
end


end
%------------- END OF CODE --------------
