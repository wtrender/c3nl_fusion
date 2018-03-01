function ix = partition(S,opt,optv,seed,g)
%% GET.PARTITION: One line description of what the function or script performs
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
%  VERSION:  0.0 CREATED: 09-Jul-2017 11:13:17
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
partition(input01,input02,input03,input04)
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
rng(seed);
if ischar(S);S = cellstr(S);end
set = unique(S);
N = numel(set);
switch opt
    case  'holdout'
        ix = randperm(N,ceil(N*(1-optv)));
        ix = ismember(S,set(ix));
    case  'gpholdout'       
        GP = unique(g);
        ix = false(size(S));
        for ii=1:numel(GP)
            s = unique(S(g==GP(ii)));
            n = numel(s);
            tix = randperm(n,ceil(n*(1-optv)));
            tix = ismember(S,s(tix));
            ix = ix|tix;
        end
        ix = logical(ix);
    otherwise
        throw(MException('Fusion:noOpt','This function has no such option'))
end


%------------- END OF CODE --------------
