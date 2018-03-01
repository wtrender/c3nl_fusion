function [valid,ext] = ext(fn,expext,err)
%% CHECK.EXT: gets a file and checks for existence and expected extension
%             if invalid throws an execption and the function is terminated  
%
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
%  VERSION:  0.0 CREATED: 24-May-2017 08:40:29
%
%% INPUTS:
%    fn - file path
%
%
%% OUTPUT:
%    valid - exists and is valid
%    ext - the extension
%
%% EXAMPLES:
%{
[valid,ext] = check.ext(fn)
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
if ~exist('err','var');err = 1;end
if exist(fn,'file')
    [~,~,ext]= fileparts(fn);
    if exist('expext','var') 
        if any(strcmpi(expext,ext))
            valid = 3;
        else 
            valid =2;
        end
    else
        valid = 1; % if a comparable extension does not exist then return 1 to indicate that the file exists
    end
else 
    valid =0;
end
if err
switch valid
    case 0;ME = MException('C3NL:noSuchFile','File %s not found',fn);
    case 2;ME = MException('C3NL:invalidExtension','File %s has a different extension from the expected %s',fn,expext);
    case {3,1};ME = 0;
end
if ME~=0                
    throw(ME)
end
end
end


%------------- END OF CODE --------------
