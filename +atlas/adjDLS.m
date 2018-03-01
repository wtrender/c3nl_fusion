function [idx,w] = adjDLS(A,s,r)
%% ATLAS.ADJDLS: One line description of what the function or script performs
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
%  VERSION:  0.1 CREATED: 25-Sep-2017 11:36:55
%
%% INPUTS:
%    A - adjecency matrix 
%    s - source point 
%    r - radius in steps  
%
%
%% OUTPUT:
%
% idx  = list of linear indexs of edges connected to source
% w = the amounts of steps from the source 

%% EXAMPLES:
%{
adjDLS(input01,input02,input03,input04)
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
ix = DLS(A,s,r);
[~,idx]=sort(ix(:,3),'descend');% sort to preserve precedence 
ix = ix(idx,:);
idx = ix(:,2);
w = r-ix(:,3);



end

function ix = DLS(A,s,r)
% this is the recursion bit
X = sparse(length(A),1);
X(s) = 1; % define source
tmp = unique([find(A(s,:)), find(A(:,s))']); % find decendents

if r>0 % if depth was not reached go deeper 
    ix = [];
    for ii=1:numel(tmp)
        ix = [ix;s,tmp(ii),r;DLS(A,tmp(ii),r-1)];
    end
else
    ix = [repmat(s,numel(tmp),1),tmp(:),repmat(r,numel(tmp),1)];
    
end

end
%------------- END OF CODE --------------
