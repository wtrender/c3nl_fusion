function output=glmnet(X,Y,varargin)
%% FIT.GLMNET: wrapps glmnet functionality
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
%  VERSION:  0.1 CREATED: 05-Jun-2017 14:51:20
%
%% TODO:
%  finish implementing all the knobs to use glmnet
%% INPUTS:
%    varargin -
%           'family' -
%             'gaussian', 'poisson'(non-negative counts),'binomial'
%             'multinomial', 'cox' or 'mgaussian'
%           'options'  - predefined glmnet set
%           'loss'     - 'deviance','class','auc','mse','mae'
%           'stability' -  here the idea is to build on group variance to
%                          estimate feature stability
%           'mc'        -  if stability estimate was chosen this allows to modify the number of monte carlo bootstraps default = 100
%           'seed'      -  assign a seed for repreducability default = 220475
%           'holdout'   -   default = 0.3
%           'outdir'    -   output directory to save results to
%           'nfolds'
%           'weights'
%           'offset'
%           'alpha'
%           'nlambda'
%           'lambda_min'
%           'lambda'
%           'standardize'
%           'intr'
%           'thresh'
%           'dfmax'
%           'pmax'
%           'exclude'
%           'penalty_factor'
%           'cl'
%           'maxit'
%           'gtype'
%           'ltype'
%           'standardize_resp'
%           'mtype'

%           'output'
%           'testset'
%
%
%% OUTPUT:
%
%% EXAMPLES:
%{
mdl=fit.glmnet(Xt(:,~ix),Yt,'mc',12)
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
disp('')
% make sure that the glmnet lib is lodede to the workspace
source = strsplit(which(['+ml/+fit/' mfilename]),'+');
addpath([source{1} '3rdParty' filesep 'glmnet_matlab']);
% check if varagin is even and if not return an error

if ~isempty(varargin)&& mod(numel(varargin),2)
    throw(MException('FUSION:Odd varargin','This function requires additional options to be supplied as pairs of parameter and value'));
else
    output = estimateNet(X,Y,varargin);
end




end

function output = estimateNet(X,Y,param)
% sets the glmnet parameters

prop = {'options','lambda','nfolds','nlambda','family','parallel','mc','loss','seed','stability','modeltype','model'};
input = param(1:2:end);
for p=prop
    ia = find(c3nl.strDetect(input,p{1}));
    if isempty(ia);ia=0;end
    switch p{1}
        case 'options';if ia;opt=param{ia*2};else; opt=glmnetSet;end
        case 'lambda';if ia;opt.lambda=param{ia*2};end
        case 'nlambda';if ia; opt.nlambda=param{ia*2};end
        case 'parallel';if ia; opt.parallel=param{ia*2};else;opt.parallel=0; end
        case 'nfolds';if ia; nfolds=param{ia*2};else; nfolds=0;end
        case 'stability';if ia; stability=param{ia*2};else; stability=0;end
        case 'mc';if ia; mc=param{ia*2};else;mc=100; end
        case 'loss';if ia; Loss=param{ia*2};end
        case 'seed';if ia; seed=param{ia*2};else;seed=220475;end
        case 'family';if ia; family=param{ia*2};else
                if numel(unique(Y))>2
                    family='multinomial';
                    if ~exist('Loss','var');Loss='class';end
                else
                    family='gaussian';
                    if ~exist('Loss','var');Loss='mse';end
                end
            end
        case 'modeltype'
            if stability;modeltype = 'stability';
            elseif nfolds ; modeltype = 'cvmdl';
            else; modeltype = 'standard';
            end
        case 'model'
            rng(seed);
            switch modeltype
                case 'standard'; output = glmnet(X,Y,family,opt);
                case 'cvmdl';    output = cvglmnet(X,Y,family,opt,Loss,nfolds);
                case  'stability'
                    cv = cvpartition(Y,'holdout',0.3);
                    ix = cv.training;
                    Xt = X(ix,:);XT = X(~ix,:);
                    Yt = Y(ix);YT = Y(~ix);
                    tmp = glmnet(Xt,Yt,family,opt); % model training data
                    acc = get.acc(YT,predict.glmnet(tmp,XT,'class'));% estimate acc on held out data
                    nc = numel(unique(Y)); % number of classes
                    [nm,nf] = size(X); % number of measures/features
                    S = false(mc,nf,opt.nlambda,nc); % stability per feature per lambda
                    oobAcc = zeros(mc,opt.nlambda);
                    oobCE = cell(mc,1);
                    mcdf = zeros(mc,opt.nlambda,nc);
                    opt.lambda = tmp.lambda; % get lambda from original estimate
                    if check.slurm
                        run.glmnet.mdl(X,Y,'stability',mc,tmp.lambda,acc,outdir); % run parallized version that saves results to output
                    else
                        for jj=1:mc
                            bix =  randsample(nm,nm,true); % random sample with replication
                            OOB=setdiff(1:nm,bix); %% The out-of-bag indices
                            % define boot data sample
                            Xb = X(bix,:);Yb = Y(bix,:);
                            Xoob = X(OOB,:);Yoob = Y(OOB,:);
                            tic;tmp1 = glmnet(Xb,Yb,family,opt);toc % model training data
                            acc1 = get.acc(Yoob,predict.glmnet(tmp1,Xoob,'class'));% estimate acc on oob
                            fs = predict.glmnet(tmp1,Xoob,'fs');
                            ce = min(get.crossEntropy(predict.glmnet(tmp1,Xoob,'response')),[],2);
                            mcdf(jj,:,:) = predict.glmnet(tmp1,Xoob,'df');
                            oobAcc(jj,:) = acc1;
                            oobCE{jj} = squeeze(ce);
                            S(jj,:,:,:) = fs;
                            if mod(jj,80);fprintf('*');else ;fprintf('-%i\n',jj);end
                        end
                        St = estimate.stability(S);
                        output.S = S;
                        output.St = St;
                        output.oobAcc = oobAcc;
                        output.oobCE = oobCE;
                        output.mcdf = mcdf;
                    end
                    
            end
    end
end
end





%------------- END OF CODE --------------
