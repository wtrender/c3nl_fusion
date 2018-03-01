function [loss,mdl,cvmdl,ix,acc]=svmecoc(X,Y,coding,cv,cp)

% Wrapper of fitcecoc with svm learners 
% created by Eyal Soreq Imperial College London c3nl
% 03/03/17
% INPUT
% =====
% X = the multivariate measurement data (m x n)
% Y = a multi class response vector (m x 1 of K>=2 classes )
% coding  = the different multiclass strategies 
%   one-versus-all = K learners
%   one-versus-one = K(K ? 1)/2 learners
%   binary complete = 2K ? 1 ? 1 learners
%   ternary complete = (3K ? 2K + 1 + 1)/2 learners
%   ordinal = K ? 1
%   dense random = ~ 10 log2K
%   sparse random = ~ 15 log2K

%                               'onevsone' or 'allpairs'
%                              'onevsall'
%                              'binarycomplete'
%                              'ternarycomplete'
%                              'ordinal'
%                              'sparserandom'
%                              'denserandom'
% cv = when it's a value it represents hold out precentege for external cross validation (default = 0.3)
% if cv is a vector it represents a cv logical partition where 1 is
% training and 0 is cv
% cp = number of k-fold partitions (default = 10)
% OUTPUT
% ======
% mdl = the overfitted model 
% cvmdl = a 10-fold cross validtaed model
% ix = a binary vector for training and testing 
% acc = perfomance structure 

if ~exist('cv','var')||isempty(cv);cv = 0.3;end
if numel(cv)==1
    cv = cvpartition(Y,'holdout',cv); % hold-out cross-validation set
    ix = cv.training;
else     
    ix = cv;
end

Xt = X(ix,:);XT = X(~ix,:);
Yt = Y(ix);YT = Y(~ix);

if ~exist('cp','var');cp = 10;end % Stratified cross-validation

L = templateSVM('Standardize',1);
mdl = fitcecoc(Xt,Yt,'Learners',L,'coding',coding);

cvmdl = crossval(mdl,'kfold',cp);

yt = predict(mdl,Xt);
yT = predict(mdl,XT);
ycvt = cvmdl.kfoldPredict;
loss = cvmdl.kfoldLoss;
if nargout>4
[acc.T.m,acc.T.a,acc.T.f1,acc.T.g,acc.T.PV,acc.T.TR]= plot.confusiongrad(YT,yT,categories(Y));
[acc.t.m,acc.t.a,acc.t.f1,acc.t.g,acc.t.PV,acc.t.TR]= plot.confusiongrad(Yt,yt,categories(Y));
[acc.cvt.m,acc.cvt.a,acc.cvt.f1,acc.cvt.g,acc.cvt.PV,acc.cvt.TR]= plot.confusiongrad(Yt,ycvt,categories(Y));
end
% cp = cvpartition(Yt,'k',10); % Stratified cross-validation
% f = @(xtr,ytr,xte,yte) confusionmat(yte,predict(mdl,xte),'order',order);
% cfMat = crossval(f,Xt,Yt,'partition',cp);
end


