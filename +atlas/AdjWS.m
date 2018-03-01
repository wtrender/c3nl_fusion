function [C,stats] = AdjWS(W,A,r,m,method)
% ATLAS.ADJWS: One line description of what the function or script performs
%
%   __           _             
%  / _|         (_)            
% | |_ _   _ ___ _  ___  _ __    
% |  _| | | / __| |/ _ \| `_ \    :- Functional and Structural 
% | | | |_| \__ \ | (_) | | | |      Integration of Neuroimages
% |_|  \__,_|___/_|\___/|_| |_|
%
%
% AUTHOR:  Eyal Soreq
%  EMAIL:  e.soreq14@imperial.ac.uk
%  AFFILIATION:  Imperial College London
%  VERSION:  0.0 CREATED: 25-Sep-2017 14:05:01
%
% INPUTS:
%    W - Vector of dissimilarity values 
%    A - Adjacency  matrix of logicals (undirected)
%    r - Radius of search 
%    m - first label integer 
%    method - maximum or minimum of peaks 
%
% OUTPUT:
%    C = a vector of labels 
%    stats = area of parcels 
%
% EXAMPLES:
%{
L = AdjWS(W,A,5,1,'max');
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
if ~exist('m','var');m=1;end
if ~exist('r','var');r=1;end
if ~exist('method','var');method='max';end
G = graph(A);
%% the simple watershed algo
switch method
    case 'min';[~,idx]=sort(W); 
    case 'max';[~,idx]=sort(-W);
    otherwise; disp('method only supports max or min')
end
N=size(idx,1);
C=zeros(size(W));
for ii=1:N
    ix = nearest(G,idx(ii),r);
    c=C(ix);
    c=c(c>0);
    if isempty(c);C(idx(ii))=m;m=m+1;
    elseif ~any(diff(c));C(idx(ii))=c(1);
    else ;C(idx(ii))=mode(c);
    end
    if ~mod(ii,10)
        if mod(ii,1000);fprintf('*');else;fprintf('\n');end
    end
end
stats = struct2array(regionprops(Lmx,{'area'}))';

end
%------------- END OF CODE --------------
