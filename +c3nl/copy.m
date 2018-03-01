function FN = copy(fn,outdir)
%% C3NL.COPY: copy fn to output dir (if outdir doesn't exit creates it)
%    __           _             
%   / _|         (_)            
%  | |_ _   _ ___ _  ___  _ __    
%  |  _| | | / __| |/ _ \| '_ \    :- Functional and Structural 
%  | | | |_| \__ \ | (_) | | | |      Integration of Neuroimages
%  |_|  \__,_|___/_|\___/|_| |_|
%
%
%% AUTHOR:  Eyal Soreq
%  EMAIL:  e.soreq14@imperial.ac.uk
%  AFFILIATION:  Imperial College London
%  VERSION:  0.0 CREATED: 24-May-2017 13:43:07
%
%% INPUTS:
%    fn - 
%    outdir - 
%
%
%% OUTPUT:
%
%% EXAMPLES:
%{
c3nl.copy(fn,outdir)
%}
%
%% DEPENDENCIES:
%
% This file is part of C^3NL Pipeline
% C^3NL Pipeline is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% C^3NL Pipeline is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with C^3NL Pipeline.If not, see <http://www.gnu.org/licenses/>.
%------------- BEGIN CODE --------------
%

new.folder([],outdir);

if iscell(fn)
    FN = {};
	for ii=1:numel(fn)
		copyfile(fn{ii},outdir);
        source = fileparts(fn{ii});
		if mod(ii,80); fprintf('*');else ;fprintf('\n');end
        FN{ii} = strep(fn{ii},source,outdir);
	end
else
    copyfile(fn,outdir);
    source = fileparts(fn);
    FN = strep(fn,source,outdir);
end

end

%------------- END OF CODE --------------
