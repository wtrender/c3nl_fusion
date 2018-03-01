function output = eco(X,Y,S,seed,ho,k)
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
RUSecoc(input01,input02,input03,input04)
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
ix = any(isnan(X));
X = X(:,~ix);
if ~exist('seed','var');seed=22041975;end
if ~exist('ho','var');ho=0.1;end
if ~exist('k','var');k=10;end

rng(seed);
ix = get.partition(S,'gpholdout',ho,seed,Y);
XR = X(~ix,:);YR = Y(~ix);
Xt = X(ix,:);Yt = Y(ix);
L =templateSVM('Standardize',1);
coding = 'onevsall';

EC = fitcecoc(Xt,Yt,'Learners',L,'coding',coding); % learn
Null = fitcecoc(Xt,Yt(randperm(size(Yt,1))),'Learners',L,'coding',coding); % estimate Null
CV = crossval(EC,'kfold', k); % estimate loss within box  

[nF1,nf1] = ml.score.F1(Yt,predict(Null,Xt));
[nF1r,nf1r] = ml.score.F1(YR,predict(Null,XR));
[F1,f1] = ml.score.F1(Yt,predict(EC,Xt));
[F1r,f1r] = ml.score.F1(YR,predict(EC,XR));
nacc = [ml.score.acc(Yt,predict(Null,Xt)),ml.score.acc(YR,predict(Null,XR))];
acc = [ml.score.acc(Yt,predict(EC,Xt)),ml.score.acc(YR,predict(EC,XR))];
results = [table({'Null'},nacc,[nF1,nF1r],nf1,nf1r,'VariableNames',{'Model','Acc','F1','f1','f1r'});
           table({'True'},acc,[F1,F1r],f1,f1r,'VariableNames',{'Model','Acc','F1','f1','f1r'})];
output = struct('mdl',EC,'loss',CV.kfoldLoss,'Performance',results);
end

    
    


%------------- END OF CODE --------------
