function output = ecoc(id,obj,perm)
%% RUN.ML.ECOC: One line description of what the function or script performs
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
%  VERSION:  0.0 CREATED: 09-Jul-2017 02:51:09
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
ecoc(input01,input02,input03,input04)
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
if isempty(id);id='setup';end
runType = find([strcmpi(id,'setup'),strcmpi(id,'clean'),isnumeric(id)]);
if ischar(obj);load(obj);else; mdl=obj;end
switch runType
    case 1;disp('Setup parallel dependencies')
        mdl = setparams(obj);
        data = ['data=' mdl.D.data filesep 'data.mat'];
        save([mdl.D.data filesep 'data.mat'],'mdl','-V7.3');
        if ~mdl.perm
            mdl.cmd = {};
            mdl.cmd{1} = struct('cmd',{{data ;'matlab -nodisplay -nodesktop -nosplash -nojvm -r "run.ml.ecoc($SLURM_ARRAY_TASK_ID,''$data'');exit"'}},...
                'jobname',mdl.name,...
                'workdir',mdl.D.tmp,...
                'output',mdl.D.output,...
                'array',sprintf('1-%i',mdl.mc),'partition','veryshort');
            mdl.cmd{2} = struct('cmd',{{data;'matlab -nodisplay -nodesktop -nosplash -r "run.ml.ecoc(''clean'',''$data'');exit"'}},...
                'jobname',['c_' mdl.name],'workdir',mdl.D.tmp,'output',mdl.D.output,'partition','veryshort');
            jid = c3nl.sbatch(mdl.cmd{1});
            mdl.cmd{2}.depend = jid;
            c3nl.sbatch(mdl.cmd{2})
        else
            mdl.cmd = {};
            for ii=1:mdl.perm
                mdl.cmd{1} = struct('cmd',{{data ;'matlab -nodisplay -nodesktop -nosplash -nojvm -r "run.ml.ecoc($SLURM_ARRAY_TASK_ID,''$data'',' num2str(ii) ');exit"'}},...
                    'jobname',mdl.name,...
                    'workdir',mdl.D.tmp,...
                    'output',mdl.D.output,...
                    'array',sprintf('1-%i',mdl.mc),'partition','veryshort');
                mdl.cmd{2} = struct('cmd',{{data;'matlab -nodisplay -nodesktop -nosplash -r "run.ml.ecoc(''clean'',''$data'');exit"'}},...
                    'jobname',['c_' mdl.name],'workdir',mdl.D.tmp,'output',mdl.D.output,'partition','veryshort');
                jid = c3nl.sbatch(mdl.cmd{1});
                mdl.cmd{2}.depend = jid;
                c3nl.sbatch(mdl.cmd{2})
            end
        end
    case 2;disp('Clean up and agragate results')
        
    case 3;fprintf('RUN job %i',id)
        mdl = setparams(mdl);
        if exist('perm','var')
            X = mdl.X(mdl.ix(:,perm),:);Y = mdl.Y(mdl.ix(:,perm));
            mdl.XT = mdl.X(~mdl.ix(:,perm),:);mdl.YT = mdl.Y(~mdl.ix(:,perm));
        else
            X = mdl.X(mdl.ix,:);Y = mdl.Y(mdl.ix);
            mdl.XT = mdl.X(~mdl.ix,:);mdl.YT = mdl.Y(~mdl.ix);            
        end
        rng(mdl.seed(id));
        [nm,~] = size(X);
        [Yb,bix] = datasample(Y,nm,1,'Replace',true);% random sample with replication
        Xb = mdl.X(bix,:);
        oobix=setdiff(1:nm,bix); %% The out-of-bag indices
        Xo = mdl.X(oobix,:);Yo = mdl.Y(oobix);%% The out-of-bag data
        EC = fitcecoc(Xb,Yb,'Learners',mdl.L,'coding',mdl.coding);
        ECNull = fitcecoc(Xb,Yb(randperm(size(Yb,1))),'Learners',mdl.L,'coding',mdl.coding);
        cvEC = crossval(EC,'kfold',mdl.k);
        yb = predict(EC,Xb);% assess overfitting
        yo = predict(EC,Xo);% assess oob
        ycvt = cvEC.kfoldPredict;% assess crossvalidation
        yon = predict(ECNull,Xo);% assess null
        output.loss = cvEC.kfoldLoss;
        output.mdl = EC;
        output.cvmdl = cvEC;
        output.perf = [get.acc(Yb,yb),get.acc(Yo,yo),get.acc(Yb,ycvt),get.acc(Yo,yon);
                       get.F1(Yb,yb),get.F1(Yo,yo),get.F1(Yb,ycvt),get.F1(Yo,yon)];
        if isfield(mdl,'XT')
            yT = predict(EC,mdl.XT);% assess Test
            output.Test = [get.acc(mdl.YT,yT),get.F1(mdl.YT,yT)];
        end
        if isfield(mdl,'XR')
            yR = predict(EC,mdl.XR);% assess Test
            output.Replicate = [get.acc(mdl.YR,yR),get.F1(mdl.YR,yR)];
        end
        if exist('perm','var')
            save(sprintf('%s%s%s_%i.%s',mdl.D.results.(sprintf('perm%03d',perm)),filesep,mdl.name,id,'mat'),'output','-V7.3');
        else
            save(sprintf('%s%s%s_%i.%s',mdl.D.results,filesep,mdl.name,id,'mat'),'output','-V7.3');
        end


    otherwise
        ME = MException('Fusion:noSuchRunID','unknown run id: %s',id);
        throw(ME)
end
end


function mdl = setparams(obj)

prop = {'name','root','results','X','Y','coding','Learners','kfold','ho','perm','mc','seed','stratify','YT','XT','YR','XR'};

if ~isstruct(obj);ia = zeros(numel(prop),1); else; ia = ismember(prop,fieldnames(obj));end
for ii=1:length(ia)
    tmp = prop{ii};
    switch tmp
        case 'name';if ia(ii);mdl.name=obj.(tmp);else; mdl.name='ecocModel';end
        case 'root';if ia(ii);mdl.root=obj.(tmp);else; mdl.root=pwd;end
        case 'results';if ia(ii);mdl.D=obj.(tmp);else; mdl.D=new.folder(mdl.name,mdl.root,{'data','results','output','tmp'});end
        case 'X';if ia(ii);mdl.X=obj.(tmp);else;throw(MException('Fusion:noData','This function requires data to run'));end
        case 'Y';if ia(ii);mdl.Y=obj.(tmp);else;throw(MException('Fusion:noResponse','This function requires a reponse vector to fit to'));end
        case 'coding';if ia(ii);mdl.coding=obj.(tmp);else; mdl.coding='onevsall';end
        case 'Learners';if ia(ii);mdl.L =obj.(tmp);else; mdl.L =templateSVM('Standardize',1);end
        case 'kfold';if ia(ii);mdl.k =obj.(tmp);else; mdl.k =3;end
        case 'perm';if ia(ii);mdl.perm=obj.(tmp);else; mdl.perm=0;end
        case 'mc';if ia(ii);mdl.mc=obj.(tmp);else; mdl.mc=1;end
        case 'ho';if ia(ii);mdl.ho=obj.(tmp);else; mdl.ho=0.25;end
        case 'seed';if ia(ii);mdl.seed=obj.(tmp);else; mdl.seed=randi([1e8,1e9],mdl.mc,1);end
        case 'stratify'
            if ia(ii)&&~isempty(obj.(tmp))
                mdl.stratify=obj.(tmp);
                if mdl.perm
                    mdl.permSeed=randi([1e8,1e9],mdl.perm,1);
                    mdl.D.results = new.folder([],mdl.D.results,arrayfun(@(x) sprintf('perm%03d',x),1:mdl.perm,'un',0));
                    ix = false(size(mdl.X,1),mdl.perm);
                    for p=1:mdl.perm
                        ix(:,p) = get.partition(mdl.stratify,'holdout',mdl.ho,mdl.permSeed(p));
                    end
                else
                    ix = get.partition(mdl.stratify,'holdout',mdl.ho,22041975);
                end
                mdl.ix=ix;
                mdl.s = mdl.stratify;
            end
            mdl.stratify = [];
        case 'YT';if ia(ii);mdl.YT=obj.(tmp);end
        case 'XT';if ia(ii);mdl.XT=obj.(tmp);end
        case 'YR';if ia(ii);mdl.YR=obj.(tmp);end
        case 'XR';if ia(ii);mdl.XR=obj.(tmp);end
    end
end


end
%------------- END OF CODE --------------
