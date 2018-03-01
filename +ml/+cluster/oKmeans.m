function [IDX, C, SUMD, D] = oKmeans(X, k)
%% C3NL.OKMEANS: Wrapper to matlab kmeans that retuns an ordered index and ordered centers 
% kmeans
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
%  VERSION:  0.0 CREATED: 09-Jun-2017 08:16:29
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
[IDX, C, SUMD, D] = c3nl.oKmeans(X, k)
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

[IDX, C, SUMD, D] = kmeans(X, k);
[C,ix]=sort(C);% reorder centers
IDX = dummyvar(IDX);
[~,IDX] = max(IDX(:,ix),[],2);
SUMD = SUMD(ix);

%------------- END OF CODE --------------
