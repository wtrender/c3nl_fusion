function [F1score,f1] = F1(Y,Yfit,class)
%% GET.F1: One line description of what the function or script performs
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
%  VERSION:  0.0 CREATED: 10-Jun-2017 23:04:57
%
%% INPUTS:
%    Y - 
%    Yfit - 
%    class - 
%
%
%% OUTPUT:
%    F1score -
%    f1 -
%% EXAMPLES:
%{

[F1score,f1] = get.F1(Y,Yfit,class)
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
if ~exist('class','var');class = categories(Yfit);end
F1score = zeros(size(Yfit,2),1);
f1 = zeros(size(Yfit,2),numel(class));
cw = sum(dummyvar(Y));
for ii=1:size(Yfit,2)
    cmat = confusionmat(Yfit(:,ii),Y,'ORDER',class)';
% True rate the proportion between the true positive and the ground truth
% aka Sensitivity Recal
    TR =  diag(cmat)./ sum(cmat,2);TR = [TR,1-TR]*100;
% Predictive Value the proportion between the true positive and the any
% value classified in class k aka Precision
    PV = diag(cmat)./sum(cmat)';PV = [PV,1-PV]*100;
    f1(ii,:) =   2*((PV(:,1).*TR(:,1))./(PV(:,1) + TR(:,1)));
    F1score(ii) = nanmean(f1(ii,:));
end
%------------- END OF CODE --------------
