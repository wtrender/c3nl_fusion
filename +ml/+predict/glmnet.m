function varargout = glmnet(mdl,X,outtype,lambda)
%% PREDICT.GLMNET: wrapper over the glmnetpredict to quicly output response vectors
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
%  VERSION:  0.0 CREATED: 05-Jun-2017 16:48:11
%
%% INPUTS:
%    mdl - a glmnet model
%    X - data that corrosponds to the model
%    type -  'link','response','coefficients','class','nonzero'
%
%    lambda - a specific lambda to get
%
%
%% OUTPUT:
%
%% EXAMPLES:
%   [yc,yp,yl]=glmnet(mdl,X)
%{

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

%% make sure that the lib is in the path

source = strsplit(which(['+ml/+predict/' mfilename]),'+');
addpath([source{1} '3rdParty' filesep 'glmnet_matlab']);

%%

resp = {'class','response','link','coefficients','nonzero','fs','dfmat','df','sign'};
if ~exist('lambda','var');lambda = mdl.lambda;end
if ~exist('outtype','var')|| isempty(outtype);outtype = {'class','response','link'};
else
    if ischar(outtype); outtype={outtype};end
end
ia = find(ismember(resp,outtype));
c = 1;
varargout = cell(nnz(ia),1);
for ii=1:numel(ia)   
    switch  resp{ia(ii)}
        case  'class';varargout{c}= glmnetPredict(mdl,X,lambda,'class');c=c+1;
        case  'response';varargout{c}= glmnetPredict(mdl,X,lambda,'response');c=c+1;
        case  'link';varargout{c}= glmnetPredict(mdl,X,lambda,'link');c=c+1;
        case  'coefficients'
            tmp = glmnetPredict(mdl,X,lambda,'coefficients');
            beta = zeros([size(tmp{1}),numel(tmp)]);
            for jj=1:numel(tmp)
                beta(:,:,jj) = tmp{jj};
            end
            varargout{c}= beta;c=c+1;
        case  'nonzero';  varargout{c}= glmnetPredict(mdl,X,lambda,'nonzero');c=c+1;
        case  'dfmat'
            tmp = predict.glmnet(mdl,X,'nonzero');
            varargout{c} = sum(cell2mat(tmp'));c=c+1;
        case  'df'
            tmp = predict.glmnet(mdl,X,'fs');
            varargout{c} = squeeze(sum(tmp));c=c+1;
        case  'sign';varargout{c} =sign(predict.glmnet(mdl,X,'coefficients'));c=c+1;
        case  'fs'
            tmp = predict.glmnet(mdl,X,'nonzero');
            fs = zeros([size(tmp{1}),numel(tmp)]);
            for jj=1:numel(tmp)
                fs(:,:,jj) = tmp{jj};
            end
            varargout{c} = fs;c=c+1;
    end
end

end


%------------- END OF CODE --------------
