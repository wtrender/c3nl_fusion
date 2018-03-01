function [idx,cost] = assignment(C,type,dir)
%% C3NL.ASSIGNMENT: Algorithm Optimal Assignment (Kuhn and Munkres)
% complexity O(n^3) https://en.wikipedia.org/wiki/Hungarian_algorithm
% adapted from http://www.mathworks.com/matlabcentral/fileexchange/20652-hungarian-algorithm-for-linear-assignment-problems--v2-3-
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
%  VERSION:  0.0 CREATED: 14-Jun-2017 02:26:54
%
%     Explenation
%     =====
%     we have a complete bipartite graph G = (S,T;E) where G is the graph
%     elment S are the source vertics and T are the target vertics E are all
%     the possible directed weigted edges.
%     our goal is to find the subgraph that produces the best matching with
%     minimum cost.
% 
%% INPUTS:
%    C - a cost matrix where rows represent set 1 and columns represent set
%        2 of a complete bipartite graph.
%    type -  the type of assignment we are after (i.e. rows or columns)
%    dir - the direction we want to go (i.e. either 'min' or 'max')
%
%
%% OUTPUT:
%    idx -  the sorting idx
%    cost -  the cost i.e. the sum(diag(M))
%
%% EXAMPLES:
%{
[idx,cost] = c3nl.assignment(M,'rows','max');
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
%% step 0 - make sure all the settings are correct
M = C;
sz = size(M); % get the size of the matrix
if (sz(1)~=sz(2)) % if the matrix is not rectangular make it so by padding inf's
    [~,ix] = max(sz);
    if ix == 1;M = padarray(M,[0 sz(1)-sz(2)],inf,'post');C=M;
    else M = padarray(M,[sz(2)-sz(1),0],inf,'post');C = M;end
end
sz = size(M);

switch dir
    case 'min';M(M(:)==inf)=max((M(M~=inf)))+std(M(M~=inf)); % make added values unatractive base on direction
    case 'max';M(M(:)==inf)=min((M(M~=inf)))-std(M(M~=inf));M = M*-1; % make padded values unatractive base on direction then inverse the matrix
    otherwise ; disp('this function only supports min and max sorting');
end

switch type
    case 'rows';M=M';disp('Assigning row''s to column''s');
    case 'cols';disp('Assigning column''s to row''s');
    otherwise;disp('only row''s or column''s are a viable selection');
end
idx = zeros(sz(1),1);
cost = 0;
n = sz(1);
%% step 1 - substract the row minimum from each row - the logic is simple removing the minimum uncovers the real overhead cost of assignment
minR = min(M,[],2); % find the minimum cost per row
% sum(minRow) will be the best cost possible in a perfect world while
% disregarding the notion that some matches are not unique
% minRed = bsxfun(@minus, M, minRow) the inner term is the matrix after we
% reduced the minRow
minC = min(bsxfun(@minus, M, minR)); % find those column where a perfect match is not present and further asses the cost of relaxing the problem
% sum(minCol) the added overhead of cost for those columns that don't have
% minimal row match
zP = M == bsxfun(@plus, minC, minR); % mark all the zeros
starZ = zeros(n,1);% get intial assignment from the Mzero stage
while any(zP(:))
    [r,c]=find(zP,1);% find one zero cost matching
    starZ(r)=c; % place it in the bucket
    zP(r,:)=false; % remove all the options of matching for both source
    zP(:,c)=false; % and target
end

while true
    %% step 2 - if M is perfect then return. otherwise goto step 3
    if all(starZ>0);break;end
    %% step 3 - Cover each column with a starred zero. If all the columns are
    %           covered then the matching is maximum if not goto step 4
    coverColumn = false(1,n);
    coverColumn(starZ(starZ>0))=true;
    coverRow = false(n,1);
    primeZ = zeros(n,1);
    [rIdx, cIdx] = find(M(~coverRow,~coverColumn)==bsxfun(@plus,minR(~coverRow),minC(~coverColumn)));% find the indexes from subgraph that isn't in starZ
    %% step 4 - Find a noncovered zero and prime it.
    while true
        cR = find(~coverRow);
        cC = find(~coverColumn);
        rIdx = cR(rIdx);
        cIdx = cC(cIdx);
        Step = 6;
        while ~isempty(cIdx)
            uZr = rIdx(1);
            uZc = cIdx(1);
            primeZ(uZr) = uZc;
            stz = starZ(uZr);
            if ~stz
                Step = 5;
                break;
            end
            coverRow(uZr) = true;
            coverColumn(stz) = false;
            z = rIdx==uZr;
            rIdx(z) = [];
            cIdx(z) = [];
            cR = find(~coverRow);
            z = M(~coverRow,stz) == minR(~coverRow) + minC(stz);
            rIdx = [rIdx(:);cR(z)];
            cIdx = [cIdx(:);stz(ones(sum(z),1))];
        end
        if Step == 6
            [minval,rIdx,cIdx]=outerplus(M(~coverRow,~coverColumn),minR(~coverRow),minC(~coverColumn));
            minC(~coverColumn) = minC(~coverColumn) + minval;
            minR(coverRow) = minR(coverRow) - minval;
        else
            break
        end
    end
    rowZ1 = find(starZ==uZc);
    starZ(uZr)=uZc;
    while rowZ1>0
        starZ(rowZ1)=0;
        uZc = primeZ(rowZ1);
        uZr = rowZ1;
        rowZ1 = find(starZ==uZc);
        starZ(uZr)=uZc;
    end
end
% Cost of assignment
idx = starZ;
switch type
    case 'rows';C_diag = diag(C(idx,:));
    case 'cols';C_diag = diag(C(:,idx));    
end
cost = sum(C_diag(C_diag~=inf))./numel(C_diag);
end

function [minval,rIdx,cIdx]=outerplus(M,x,y)
ny=size(M,2);
minval=inf;
for c=1:ny
    M(:,c)=M(:,c)-(x+y(c));
    minval = min(minval,min(M(:,c)));
end
[rIdx,cIdx]=find(M==minval);
end
%------------- END OF CODE --------------
