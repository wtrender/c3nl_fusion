function output = MCglmnet(X,Y,S,seed,lf,ls)
%% RUN.ML.MCglmnet: One line description of what the function or script performs
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
%  VERSION:  0.0 CREATED: 29-Nov-2017 10:44:37
%
%% INPUTS:
%    X - 
%    Y - 
%    seed - 
%
%
%% OUTPUT:
%
%% EXAMPLES:
%{
ml.fit.RUSglmnet(input01,input02,input03,input04)
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

% for each class random sample with replication to reaching the smallest
% group number of measures 
if ~exist('lf','var');lf=1;end
if ~exist('ls','var');ls=1;end

ixn = any(isnan(X),2);
X = X(~ixn,:);
Y = Y(~ixn);
S = S(~ixn);
%L = L(~ixn);
warning off
if ~exist('seed','var');seed=22041975;end
rng(seed);
ix = get.partition(S,'holdout',0.3,seed);

XR = X(~ix,:);YR = Y(~ix);
Xt = X(ix,:);Yt = Y(ix);
f = @(t,mx,mn,a) exp((t./numel(t)).*log(a*mn/mx))./mn;
lambda = f(1:100,size(Xt,2),numel(unique(Yt)),lf);
try
    cv = ml.fit.glmnet(Xt,Yt,'family','multinomial','nfolds',10,'loss','class','lambda',lambda);
catch 
    cv.glmnet_fit = ml.fit.glmnet(Xt,Yt,'family','multinomial','lambda',lambda);
    ls = 3;
    cv.cvm = zeros(size(lambda));
    cv.nzero = cv.glmnet_fit.df;
end
nc = numel(unique(Y));

ce = cv.cvm;
df = cv.nzero;
yp = squeeze(mean(min(ml.score.crossEntropy(ml.predict.glmnet(cv.glmnet_fit,Xt,'response')),[],2)));
switch ls
    case 1; idx = find(lambda<=cv.lambda_1se,1,'first');
    case 2; idx = find(lambda<=cv.lambda_min,1,'first');
    case 3; idx = plot.scree(yp);
    case 4; idx = plot.scree(ce);
end
%idx = plot.scree(yp);
%idx1 = find(ce<(ce(idx)-std(ce)),1,'last');
%[~,idx2] = min(ce);
mdl = ml.fit.glmnet(Xt,Yt,'family','multinomial','lambda',lambda);
null = ml.fit.glmnet(Xt,Yt(randperm(size(Yt,1))),'family','multinomial','lambda',lambda);

[nF1r,nf1r] = ml.score.F1(YR,ml.predict.glmnet(null,XR,[],mdl.lambda(idx)));
[F1r,f1r] = ml.score.F1(YR,ml.predict.glmnet(mdl,XR,[],lambda(idx)));
nacc = [ml.score.acc(YR,ml.predict.glmnet(null,XR,[],lambda(idx)))];
acc = [ml.score.acc(YR,ml.predict.glmnet(mdl,XR,[],lambda(idx)))];
Beta = zeros(size(X,2),nc);
for ii=1:nc
Beta(:,ii) = mdl.beta{ii}(:,idx);
end
results = [table({'Null'},nacc,nF1r,nf1r,'VariableNames',{'Model','Acc','F1','f1r'});
           table({'True'},acc,F1r,f1r,'VariableNames',{'Model','Acc','F1','f1r'})];
output = struct('mdl',mdl,'loss',ce(idx),'Performance',results,'W',Beta,'df',nnz(any(Beta,2)));
end

    
    


%------------- END OF CODE --------------