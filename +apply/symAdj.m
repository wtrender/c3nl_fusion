function Adj = symAdj(Adj,type)
%% C3NL.SYMADJ: symmetrize a square matrix using different ways
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
%  VERSION:  0.0 CREATED: 07-Jun-2017 16:31:54
%
%% INPUTS:
%    Adj - 
%    type - 
%
%% OUTPUT:
%
%% EXAMPLES:
%{
c3nl.symAdj(Adj)
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
A =     Adj;
AA =    fliplr(rot90(Adj,-1));
ix_l = tril(ones(size(Adj)),-1)>0;
ix_u = triu(ones(size(Adj)),1)>0;
switch type
    case 'lower'; Adj = ix_l.*A+ix_u.*AA;
    case 'upper'; Adj = ix_u.*A+ix_l.*AA;
    case 'mean';  Adj = mean(cat(3,c3nl.symAdj(Adj,'lower'),c3nl.symAdj(Adj,'upper')),3);
    otherwise; disp('type is not recognized');return;
end
end
%------------- END OF CODE --------------
