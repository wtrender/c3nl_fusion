function [l,t] = reLabel(L,offset,T)
%% ATLAS.RELABEL: One line description of what the function or script performs
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
%  VERSION:  0.0 CREATED: 08-Aug-2017 14:24:47
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
reLabel(input01,input02,input03,input04)
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
if ~exist('offset','var');offset=0;end
idL = unique(L(~isnan(L))); % you atlas id's
idL(idL==0)=[];% ignore zero's
nL = numel(idL);
if max(idL)>=nL
    LUT = zeros(2^16,1,'uint16');
    LUT(idL+1) = offset+1:offset+nL;
    l = intlut(uint16(L),LUT);
else 
    l=L;
end

l(isnan(L))=nan;
if exist('T','var') && nargout>1
    t = T;
    t.oldROIid = T.ROIid;
    t.ROIid = (offset+1:offset+nL)';
else
    t = table((offset+1:offset+nL)',idL);
    t.Properties.VariableNames = {'ROIid','oldROIid'};
end
%------------- END OF CODE --------------
