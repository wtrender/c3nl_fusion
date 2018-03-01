function [D,Group] = project(root,projectname,folders)
%% NEW.PROJECT: either creates a new project structure or loads the current 
%               group and logs the entry in a text log file 
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
%  VERSION:  0.0 CREATED: 23-May-2017 11:42:55
%
%% INPUTS:
%    inDir - home directory where the project will be placed in
%    projectname - the project name   
%    folders - folders to create in the project folder default are {'tmp','output','LOG','GROUP'}
%
%
%% OUTPUT:
%    D - a struct of strings pointing at the full path's of the subfolders 
%    Group - struct containg group info about the experiment 
%
%% EXAMPLES:
%{
[D,group] = project(inDir,folders)
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
if ~exist('folders','var')
    folders = {'tmp','output','LOG','GROUP'};
end
if ~ismember('GROUP',folders);folders{end+1} = 'GROUP';end
D = new.folder([],root,folders);

fid = new.file('log.txt',D.GROUP,'a+');
fprintf(fid,'%s\n',datestr(datetime('now')));
if ~(exist([D.GROUP,filesep,'Group.mat'],'file'))
    disp('Constructing group log structure')
    Group = [];
    Group.D = D;
    Group.name =projectname;
    Group.log = [D.GROUP,filesep,'log.txt'];
    fn = fieldnames(D);
    for ii=1:numel(D);fprintf(fid,'%s\n',D.(fn{ii}));end
    save ([D.GROUP,filesep,'Group.mat'], 'Group');
else 
    load([D.GROUP,filesep,'Group.mat']);
end
end
    
%------------- END OF CODE --------------
