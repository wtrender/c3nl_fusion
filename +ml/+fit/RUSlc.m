function output = RUSlc(X,Y,S,seed)
%% RUN.ML.RUSECOC: One line description of what the function or script performs
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
ml.fit.RUSecoc(input01,input02,input03,input04)
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

ixn = any(isnan(X),2);
X = X(~ixn,:);
Y = Y(~ixn);
S = S(~ixn);
%L = L(~ixn);

if ~exist('seed','var');seed=22041975;end
rng(seed);
ix = get.partition(S,'holdout',0.3,seed);
XR = X(~ix,:);YR = Y(~ix);
Xt = X(ix,:);Yt = Y(ix);
k=7;
nu = min(sum(dummyvar(Yt)));
nc = unique(Y);
bix = [];% generate balanced box set for training without replacment 
oix = [];% generate imbalanced out-of-box set for validation without replacment  
for ii=1:numel(nc)
    ix = find(Yt==nc(ii));
    [~,tix] = datasample(Yt(ix),nu,1,'Replace',false);
    bix = [bix;ix(tix)];
    dix=setdiff(1:numel(ix),tix); %% The out-of-bag indices
    oix = [oix;ix(dix)];
end
Xb = Xt(bix,:); % sampled data
Yb = Yt(bix); % sampled labels
Lambda = logspace(-6,-0.5,11);
cv = fitclinear(Xb,Yb,'ObservationsIn','rows','KFold',5,...
    'Learner','logistic','Solver','sparsa','Regularization','lasso',...
    'Lambda',Lambda,'GradientTolerance',1e-8,'ScoreTransform','logit');

ce = kfoldLoss(cv);
idx = [find(diff(ce),1,'first')-1,find(diff(ce),1,'last')+1];
ul = log10(Lambda(idx));
Lambda = logspace(ul(1),ul(2),100);
cv = fitclinear(Xb,Yb,'ObservationsIn','rows','KFold',5,...
    'Learner','logistic','Solver','sparsa','Regularization','lasso',...
    'Lambda',Lambda,'GradientTolerance',1e-8,'ScoreTransform','logit');
ce = kfoldLoss(cv);
df=  sum(cv.Trained{1}.Beta~=0);
idx = find(ce<(min(ce)+std(c163e)/10),1,'first');
%idx = plot.scree(1./(df+1).*ce,1);
%idx1 = find(ce<(ce(idx)-std(ce)),1,'last');
%[~,idx2] = min(ce);
mdl = fitclinear(Xb,Yb,'ObservationsIn','rows',...
    'Learner','logistic','Solver','sparsa','Regularization','lasso',...
    'Lambda',Lambda(idx),'GradientTolerance',1e-8,'ScoreTransform','logit');

null = fitclinear(Xb,Yb(randperm(size(Yb,1))),'ObservationsIn','rows',...
    'Learner','logistic','Solver','sparsa','Regularization','lasso',...
    'Lambda',Lambda(idx),'GradientTolerance',1e-8,'ScoreTransform','logit');

[nF1r,nf1r] = ml.score.F1(YR,predict(null,XR));
[F1r,f1r] = ml.score.F1(YR,predict(mdl,XR));
nacc = [ml.score.acc(YR,predict(null,XR))];
acc = [ml.score.acc(YR,predict(mdl,XR))];
results = [table({'Null'},nacc,nF1r,nf1r,'VariableNames',{'Model','Acc','F1','f1r'});
           table({'True'},acc,F1r,f1r,'VariableNames',{'Model','Acc','F1','f1r'})];
output = struct('mdl',mdl,'null',null,'loss',ce(idx),'Performance',results,'W',mdl.Beta,'df',df(idx));
end

    
    


%------------- END OF CODE --------------
