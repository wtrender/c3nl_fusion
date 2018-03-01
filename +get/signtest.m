function tbl = signtest(x,g)
%% GET.RANKSUM: One line description of what the function or script performs
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
%  VERSION:  0.0 CREATED: 23-Dec-2017 11:26:40
%
%% INPUTS:
%    x - 
%    y - 
%
%% OUTPUT:
%
%% EXAMPLES:
%{
[p,h,tbl] = get.ranksum(x,y);
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
tbl = table();
gp = unique(g);
for ii=1:numel(gp)   
[p,h,stats] = signtest(x(g==gp(ii)),0,'method','approximate');
N = sum(g==gp(ii));
tbl = [tbl;table({char(gp(ii))},N,stats.zval,p,stats.zval/sqrt(N),'VariableNames',{'Group','N','Z','p_value','Effect_Size'})];
end

end