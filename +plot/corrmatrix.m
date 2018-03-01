function corrmatrix(X,cmap,lim)
%% PLOT.CORRMATRIX: One line description of what the function or script performs
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
%  VERSION:  0.0 CREATED: 01-Mar-2018 09:56:11
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
corrmatrix(input01,input02,input03,input04)
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
cmap =jet
disp('')
d = size(X);
ix = find(fliplr(eye(d(2))));
[x,y,w,h] = get.fancyGrid(ones(d(2),1)*d(2),0.01,0.1,'same');
clf
for ii=1:numel(x)
   ax.(sprintf('p%i',ii)) = axes('Position',  [x(ii) y(ii) w h]);%axis('off');
   id = ismember(ix,ii);
   if any(id)
    [f,xi] = ksdensity(X(:,find(id)));%,'Support',lim
    patch(xi,f,cmap(ii,:),'FaceAlpha',0.5,'parent',ax.(['p' num2str(ii)]));
   else
    [I,J] = ind2sub([d(2),d(2)],ii);
    
   end
end

%------------- END OF CODE --------------
