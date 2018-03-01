function strippedT1 = skullStrip(T1,c1,c2,c3)
%% RUN.C3NL.SKULLSTRIP: One line description of what the function or script performs
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
%  VERSION:  0.0 CREATED: 15-Jan-2018 13:28:27
%
%% INPUTS:
%    T1 - T1 volume
%    c1 - gray mater mask
%    c2 - white mater mask
%    c3 - CSF mask
%
%
%% OUTPUT:
%
%% EXAMPLES:
%{
skullStrip(input01,input02,input03,input04)
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
pth = fileparts(T1);
[T1,G] = load.vol(T1);
GM = load.vol(c1);
WM = load.vol(c2);
CSF = load.vol(c3);
mask = imfill(sum(cat(4,GM,WM,CSF)>0.85,4),'holes');
strippedT1 = sum(cat(4,GM,WM,CSF)>0.75,4).*T1;
save.vol(strippedT1,G,[pth,filesep,'strippedT1.nii'],'float32');





%------------- END OF CODE --------------
