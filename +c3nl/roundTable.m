function D = roundTable(T,f)
%% C3NL.ROUNDTABLE: One line description of what the function or script performs
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
%  VERSION:  0.0 CREATED: 26-Oct-2017 00:05:22
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
roundTable(input01,input02,input03,input04)
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

disp('')
d = size(T);
D = table2cell(T);
ix = find(cellfun(@isempty,D));
if ~isempty(ix)
    for ii=1:numel(ix)
        D{ix(ii)}=nan;
    end
end
for ii=1:d(1)
    for jj=1:d(2)
        tp = find([isnumeric(D{ii,jj}),ischar(D{ii,jj}),iscell(D{ii,jj}),1],1,'first');
        switch tp
            case 1;D{ii,jj} = ceil(D{ii,jj}*10^f)/10^f;
            case 2;D{ii,jj};
            case 3;D{ii,jj}=D{ii,jj}{1};
            otherwise;D{ii,jj} = T{ii,jj};
        end
    end
end



D = cell2table(D);
if ~isempty(T.Properties.RowNames)
    D = [T.Properties.RowNames,D];
    D.Properties.VariableNames(2:end) = T.Properties.VariableNames;
else 
    D.Properties.VariableNames = T.Properties.VariableNames;
end

if strcmp(D.Properties.VariableNames{1},'Var1');D.Properties.VariableNames{1}='Terms';end
end



%------------- END OF CODE --------------
