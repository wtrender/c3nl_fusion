function output = mdl(id,obj)
%% RUN.GLMNET.MDL: parallized wrapper to glmnet modeling
%                  1. creates a tmp folder
%                  2. save the data to model with random seeds
%                  3. calls two slurm jobs
%                      a. fits a model based on slurm id and saves the data
%                      b. integrator loads all data, deletes files and saves one structure
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
%  VERSION:  0.0 CREATED: 06-Jun-2017 11:15:54
%
%% INPUTS:
%    id = a number reffeing to the current mdl to fit or an empty arry to intilize parallelization
%    obj - either a path to a data structure or a data structure
%    containing:
%        results - the folder where the data will be outputed
%        X - multivariate data to model
%        Y - ground truth to model
%        type  - 'standard','cvmdl','stability'
%        mc - number of bootstraps
%        ho - held out data to test model
%        stratify - a vector to form the partitioning should be identical to the number of events (only applicable if ho exists)
%        opt - glmnet option structure
%        nfolds - default is not to use cvglmnet if this is set by a number and additional nfold models will be fit in each step ;
%        seed - either a vector or a single number to define randomness
%
%
%% OUTPUT:
%
%% EXAMPLES:
%{
opt = glmnetSet;
run.glmnet.mdl([],struct('name','dense','X',study.X(:,ix),'Y',study.Y,'type','stability','mc',50));
run.glmnet.mdl('setup',obj)
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

source = strsplit(which(['run.glmnet.' mfilename]),'+');
addpath([source{1} '3rdParty' filesep 'glmnet_matlab']);
%%

if isempty(id);id='setup';end
runType = find([strcmpi(id,'setup'),strcmpi(id,'clean'),isnumeric(id)]);
% if the data is in a file load it
if ischar(obj);load(obj);end
% setup job based on type
switch runType
    case 1;disp('Setup parallel dependencies')
        % here we setup two bash files one a job array that sets up the
        % unique modeling. The other is a cleanup step that loads the array results
        % agregates them together and removes the unecessery files.
        mdl = setparams(obj);
        data = ['data=' mdl.D.data filesep 'data.mat'];
        mdl.cmd = {};
        mdl.cmd{1} = struct('cmd',{{data ;'matlab -nodisplay -nodesktop -nosplash -r "run.glmnet.mdl($SLURM_ARRAY_TASK_ID,''$data'');exit"'}},...
            'jobname',mdl.name,...
            'workdir',mdl.D.tmp,...
            'output',mdl.D.output,...
            'array',sprintf('1-%i',mdl.mc),'partition','veryshort');
        mdl.cmd{2} = struct('cmd',{{data;'matlab -nodisplay -nodesktop -nosplash -r "run.glmnet.mdl(''clean'',''$data'');exit"'}},...
            'jobname',['c_' mdl.name],'workdir',mdl.D.tmp,'output',mdl.D.output,'partition','veryshort');
        save([mdl.D.data filesep 'data.mat'],'mdl','-V7.3');
        jid = c3nl.sbatch(mdl.cmd{1});
        mdl.cmd{2}.depend = jid;
        c3nl.sbatch(mdl.cmd{2})
        
    case 2;disp('Clean up and agragate results')
        fn = c3nl.select(struct('pth',mdl.D.results,'name','*.mat'));
        load(fn{1});% load one to preallocate
        nm = numel(fn);
        [nf,nl,nc] = size(output.fs);
        acc = zeros(nl,nm);
        nul = zeros(nl,nm);
        ce = zeros(nl,nm);
        df = zeros(nl,nc,nm);
        f1 = zeros(nl,nc,nm);
        F1 = zeros(nl,nm);
        fs = false(nf,nl,nc,nm);
        for ii=1:numel(fn)
            load(fn{ii})
            acc(:,ii) = output.acc_o;
            F1(:,ii) = output.F1score;
            f1(:,:,ii) = output.f1;
            nul(:,ii) = output.acc_n;
            ce(:,ii) = output.ce_o;
            df(:,:,ii) = output.df;
            fs(:,:,:,ii) = output.fs;
            if mod(ii,80);fprintf('*');else; fprintf('-%i\n',ii);end
        end
        S = estimate.stability(fs);
        S.ce = ce;
        S.df = df;
        S.acc = acc;
        if isfield(mdl,'XT')
            [~,idx]=max(sum(S.fs));%find out where is the most stable largest set of features
            S.sfs = S.fs(:,idx);
            S.lambda = mdl.gfunc(1:100,1/sqrt(nnz(S.sfs)),1/nc,nc);% estimate a new lambda range
            tic;tmp1 = glmnet(mdl.X,mdl.Y,mdl.family,mdl.opt);toc % model all features
            mdl.opt.lambda = S.lambda;
            tic;tmp = glmnet(mdl.X(:,S.sfs),mdl.Y,mdl.family,mdl.opt);toc % model stable features only
            S.acc_sT = get.acc(mdl.YT,predict.glmnet(tmp,mdl.XT(:,S.sfs)));
            S.acc_fT = get.acc(mdl.YT,predict.glmnet(tmp1,mdl.XT));
            S.full_mdl = tmp1;
            S.stable_mdl = tmp;
        end
        c3nl.select(struct('pth',mdl.D.results,'name','*.mat','delete',1));
        save(sprintf('%s%s%s_%s.%s',mdl.D.results,filesep,mdl.name,id,'mat'),'S','-V7.3');

    case 3;fprintf('RUN job %i',id)
        mdl = setparams(mdl);
        rng(mdl.seed(id));
        [nm,~] = size(mdl.X);
        [Yb,bix] = datasample(mdl.Y,nm,1,'Replace',true);% random sample with replication
        Xb = mdl.X(bix,:);
        oobix=setdiff(1:nm,bix); %% The out-of-bag indices
        Xo = mdl.X(oobix,:);Yo = mdl.Y(oobix);%% The out-of-bag data
        tic;tmp = glmnet(Xb,Yb,mdl.family,mdl.opt);toc % model boot data
        tic;tmp1 = glmnet(Xb,Yb(randperm(size(Yb,1))),mdl.family,mdl.opt);toc % model null
        Yc = predict.glmnet(tmp,Xo);
        output.acc_o = get.acc(Yo,Yc); % assess misclassification accuracy on oob
        [output.F1score,output.f1] = get.F1(Yo,Yc);
        Yc = predict.glmnet(tmp1,Xo);
        output.acc_n = get.acc(Yo,Yc); % assess chance
        [output.nF1score,output.nf1] = get.F1(Yo,Yc);
        output.ce_o = mean(squeeze(min(get.crossEntropy(predict.glmnet(tmp,Xo,'response')),[],2)));% assess crossEntropy
        output.ce_n = mean(squeeze(min(get.crossEntropy(predict.glmnet(tmp1,Xo,'response')),[],2)));% assess crossEntropy
        output.df = predict.glmnet(tmp,Xo,'df');
        output.fs = predict.glmnet(tmp,Xo,'fs');
        output.mdl = tmp;
        save(sprintf('%s%s%s_%i.%s',mdl.D.results,filesep,mdl.name,id,'mat'),'output','-V7.3');
    otherwise
        ME = MException('Fusion:noSuchRunID','unknown run id: %s',id);
        throw(ME)
end


end



function mdl = setparams(obj)
prop = {'name','root','results','X','Y','family','gridfunc','nlambda','lambda','ho','mc','seed','stratify','YT','XT','type','opt','alpha','nfolds','obj'};
if ~isstruct(obj);ia = zeros(numel(prop),1); else; ia = ismember(prop,fieldnames(obj));end
for ii=1:length(ia)
    tmp = prop{ii};
    switch tmp
        case 'opt';if ia(ii);mdl.opt=obj.(tmp);else; mdl.opt=glmnetSet;end
        case 'root';if ia(ii);mdl.root=obj.(tmp);else; mdl.root=pwd;end
        case 'alpha';if ia(ii);mdl.opt.alpha=obj.(tmp);else; mdl.opt.alpha=0.95;end % forcing some elastic net
        case 'X';if ia(ii);mdl.X=obj.(tmp);else;throw(MException('Fusion:noData','This function requires data to run'));end
        case 'Y';if ia(ii);mdl.Y=obj.(tmp);else;throw(MException('Fusion:noResponse','This function requires a reponse vector to fit to'));end
        case 'YT';if ia(ii);mdl.YT=obj.(tmp);end
        case 'XT';if ia(ii);mdl.XT=obj.(tmp);end
        case 'type';if ia(ii);mdl.type=obj.(tmp);else; mdl.type='standard';end
        case 'mc';if ia(ii);mdl.mc=obj.(tmp);else; mdl.mc=1;end
        case 'ho';if ia(ii);mdl.ho=obj.(tmp);else; mdl.ho=0.2;end
        case 'name';if ia(ii);mdl.name=obj.(tmp);else; mdl.name='glmnetmodel';end
        case 'family';if ia(ii);mdl.family=obj.(tmp);else; mdl.family='multinomial';end
        case 'gridfunc'; if ia(ii);mdl.gfunc = obj.(tmp);else;mdl.gfunc = @(t,mx,mn,tx) mn.*exp((t./tx).*log(mx/mn));end
        case 'nlambda'; if ia(ii);mdl.opt.nlambda = obj.(tmp);else;mdl.opt.nlambda=100 ;end
        case 'lambda'
            if ia(ii);mdl.opt.lambda = obj.(tmp);
            else;mdl.opt.lambda=mdl.gfunc(1:mdl.opt.nlambda,1/sqrt(size(mdl.X,2)),1/numel(unique(mdl.Y)),numel(unique(mdl.Y))); ;
            end
        case 'stratify'
            if ia(ii)&&~isempty(obj.(tmp))
                mdl.stratify=obj.(tmp);
                rng(22041975);
                cv = cvpartition(mdl.stratify,'holdout',mdl.ho);
                ix = cv.training;
                mdl.XT = mdl.X(~ix,:);mdl.YT = mdl.Y(~ix);
                mdl.X = mdl.X(ix,:);mdl.Y = mdl.Y(ix);
                mdl.cv=cv;
                mdl.ix=ix;
            end
            mdl.stratify = [];
        case 'nfolds';if ia(ii);mdl.nfolds=obj.(tmp);else; mdl.nfolds=1;end
        case 'results';if ia(ii);mdl.D=obj.(tmp);else; mdl.D=new.folder(mdl.name,mdl.root,{'data','results','output','tmp'});end
        case 'seed';if ia(ii);mdl.seed=obj.(tmp);else; mdl.seed=randi([1e8,1e9],mdl.mc,1);end
    end
end
end

%------------- END OF CODE --------------
