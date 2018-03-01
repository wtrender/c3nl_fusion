function [G,L,T] = labeltree(X,conn,gp)
   
n = size(X,1);
nleaves = n+1;
T = nan(n+1,1);
L = zeros(n+1,1);
G = zeros(n+1,1);
% Each cut potentially yields an additional cluster
todo = true(n,1);

% Define cluster numbers for each side of each non-leaf node
clustlist = reshape(1:2*n,n,2);
level=0;
% Propagate cluster numbers down the tree
while(any(todo))
   % Work on rows that are now split but not yet processed
   rows = find(todo & ~conn);
   if isempty(rows), break; end
   level = level+1;
   for j=1:2    % 1=left, 2=right
      children = X(rows,j);
   
      % Assign cluster number to child leaf node
      leaf = (children <= nleaves);
      if any(leaf)
         T(children(leaf)) = clustlist(rows(leaf),j);
         L(children(leaf)) = level;
         G(children(leaf)) = gp(rows(leaf));
      end
   
      % Also assign it to both children of any joined child non-leaf nodes
      joint = ~leaf;
      joint(joint) = conn(children(joint)-nleaves);
      if any(joint)
         clustnum = clustlist(rows(joint),j);
         childnum = children(joint) - nleaves;
         clustlist(childnum,1) = clustnum;
         clustlist(childnum,2) = clustnum;
         conn(childnum) = 0;
         L(childnum) = level;
      end
   end

   % Mark these rows as done  
   todo(rows) = 0;
end
end



