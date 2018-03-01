function L = dbscan(Y,epsilon,minVoxels,r,lm,lambda,theta)
%% ATLAS.DBSCAN: One line description of what the function or script performs
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
%  VERSION:  0.0 CREATED: 10-Oct-2017 09:02:25
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
dbscan(input01,input02,input03,input04)
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
tic;
Y = double(Y);
if ~exist('r','var');r=1;end
if ~exist('lambda','var');lambda=1;end
if ~exist('theta','var');theta=1;end
if ~exist('lm','var');lm='euc';end
dim = size(Y);
nd = numel(dim);
Y = padarray(Y,ones(numel(dim),1)*(r+1),-inf); % pad image to account for bounderies
Pdim =size(Y);
[d,w] = c3nl.neighb(Pdim,r,lm); % get connected components linear index
ix = find(isnan(Y));
Y(ix) = -inf;
[~,idx] = sort(-Y(:)); % force minima
N=size(Y(:),1);
%idx=1:N;
L = zeros(Pdim);
T = zeros(Pdim);
visited=false(N,1);
noise=false(N,1);
visited(idx(Y(idx)==-inf)) = true;
m=1;
for ii=1:N
    if ~visited(idx(ii))
        visited(idx(ii)) = true;
        idx1=idx(ii);
        try;Neighbors=RegionQuery(idx1);
        catch;keyboard;end    
        if numel(Neighbors)<minVoxels
            % X(i,:) is NOISE
            noise(idx1)=true;
        else
            m=m+1;
            ExpandCluster(idx1,Neighbors,m);
        end
    end
end

    function ExpandCluster(idx,Neighbors,m)
        L(idx) = m;
        k = 1;
        while true
            j = Neighbors(k);
            
            if ~visited(j)
                visited(j)=true;
                Neighbors2=RegionQuery(idx);
                if numel(Neighbors2)>=minVoxels
                    Neighbors=[Neighbors Neighbors2];   %#ok
                end
            end
            if L(j)==0
                L(j)=m;
            end
            
            k = k + 1;
            if k > numel(Neighbors)
                break;
            end
        end
    end

    function Neighbors=RegionQuery(idx)
        idx2=idx+d;
        tmp = [(1-c3nl.scale(Y([idx;idx2])))*lambda,theta*c3nl.scale([0;w])];
        D = pdist2(tmp,tmp);
        Neighbors=idx2(D(1,2:end)<=epsilon);
    end


c = repmat({':'},nd-1,1);
for ii=1:nd % crop image back to normal
    L = shiftdim(L,1);
    L = L(1+(r+1):end-(r+1),c{:});
end
toc

end
%------------- END OF CODE --------------
