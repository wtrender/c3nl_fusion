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
            mdl.cmd{1} = struct('cmd',{{data ;'matlab -nodisplay -nodesktop -nosplash -r "run.ml.ecoc($SLURM_ARRAY_TASK_ID,''$data'');exit"'}},...
                'jobname',mdl.name,...
                'workdir',mdl.D.tmp,...
                'output',mdl.D.output,...
                'array',sprintf('1-%i',mdl.mc),'partition','verylong');
            mdl.cmd{2} = struct('cmd',{{data;'matlab -nodisplay -nodesktop -nosplash -r "run.ml.ecoc(''clean'',''$data'');exit"'}},...
                'jobname',['c_' mdl.name],'workdir',mdl.D.tmp,'output',mdl.D.output,'partition','veryshort');
            jid = c3nl.sbatch(mdl.cmd{1});
            mdl.cmd{2}.depend = jid;
            c3nl.sbatch(mdl.cmd{2})
        else
            mdl.cmd = {};
            mdl.cmd{1} = struct('cmd',{{data ;['matlab -nodisplay -nodesktop -nosplash -r "run.ml.ecoc($SLURM_ARRAY_TASK_ID,''$data'',$SLURM_ARRAY_TASK_ID);exit"']}},...
                'jobname',mdl.name,...
                'workdir',mdl.D.tmp,...
                'output',mdl.D.output,...
                'array',sprintf('1-%i',mdl.perm),'partition','verylong','mem','15GB');
            mdl.cmd{2} = struct('cmd',{{data;'matlab -nodisplay -nodesktop -nosplash -r "run.ml.ecoc(''clean'',''$data'');exit"'}},...
                'jobname',['c_' mdl.name],'workdir',mdl.D.tmp,'output',mdl.D.output,'partition','verylong','mem','5GB');
            jid = c3nl.sbatch(mdl.cmd{1});
            mdl.cmd{2}.depend = jid;
            c3nl.sbatch(mdl.cmd{2})
            
        end
    case 2;disp('Clean up and agragate results')
        fn = c3nl.select(struct('pth',mdl.D.results.root,'name','*.mat'));
        header = {'Null','boot','oob','kFoldcv','heldoutTest'};
        F1 = table();
        Acc = table();
        for ii=1:numel(fn)
            load(fn{ii})
            tmp = array2table(output.F1,'VariableNames',header);
            tmp.Replicate = output.Replicate(:,2);
            F1 = [F1;tmp];
            tmp = array2table(output.Acc,'VariableNames',header);
            tmp.Replicate = output.Replicate(:,1);
            Acc = [Acc;tmp];
        end
        out = struct('F1',F1,'Acc',Acc);
        if ~isempty(mdl.subgroup)
            subheader = cellstr(unique(mdl.subgroup.S));
            nh = cellfun(@(x) matlab.lang.makeValidName(x),cellstr(categorical(subheader).*categorical(header)),'un',0)';
            SGF1 = zeros(numel(output.sgF1)*numel(fn),numel(nh));
            SGAcc = zeros(numel(output.sgF1)*numel(fn),numel(nh));
            if isfield(mdl.subgroup,'SR')
                subheaderR = cellstr(unique(mdl.subgroup.S));
                nhR = cellfun(@(x) matlab.lang.makeValidName(x),cellstr(categorical(subheaderR).*categorical({'Acc','F1'})),'un',0)';
                SGR = zeros(numel(output.sgr)*numel(fn),numel(nhR));
            end
            for ii=1:numel(fn)
                load(fn{ii})
                for jj=1:numel(output.sgF1)
                    idx = (ii-1)*numel(output.sgF1)+jj;
                    SGF1(idx,:) = output.sgF1{jj}(:)';
                    SGAcc(idx,:) = output.sgAcc{jj}(:)';
                    if isfield(mdl.subgroup,'SR')
                        SGR(idx,:) = output.sgr{jj}(:)';
                    end
                end
            end
            out.SGF1 = array2table(SGF1,'VariableNames',nh);
            out.SGAcc = array2table(SGAcc,'VariableNames',nh);
            if isfield(mdl.subgroup,'SR')
                out.SGR = array2table(SGR,'VariableNames',nhR(:));
            end
        end
        
        if isfield(mdl,'XT')
            tmp = [];
            nhR = cellfun(@(x) matlab.lang.makeValidName(x),cellstr(unique([mdl.subgroup.S;mdl.subgroup.ST]).*categorical({'Acc','F1'})),'un',0)';
            for ii=1:numel(fn)
                load(fn{ii});
                tmp = [tmp;[cell2mat(output.perf_t(:,1)),cell2mat(cellfun(@(x) [x(1,:),x(2,:)],output.perf_t(:,2),'un',0))]];
            end
            out.hierarchy = array2table(tmp,'VariableNames',nhR);
        end
        
        %system(['rm -rf ' mdl.D.results.root filesep '*']);
        save(sprintf('%s%s%s_%s.%s',mdl.D.results.root,filesep,mdl.name,'cleaned','mat'),'out','-V7.3');
        
        
    case 3;fprintf('RUN job %i',id)
        mdl = setparams(mdl);
        mdl.permid = perm;
        [X,XT,XR,Y,YT,YR,S,SR,ST,L,k,partmp,Replicate,seed,coding,nm,loss,F1,Acc,gF1,gAcc,gRF1,gRAcc] = extractparapms(mdl);
        tic
        C = unique(Y);
        fs = zeros(size(XT,2),numel(unique(Y)),mdl.mc);% features models mc
        if ~isempty(mdl.subgroup)  
            fsg = zeros(size(XT,2),numel(unique(Y)),numel(unique(S))+1,mdl.mc);% features models subgroups mc
        end
        for idx = 1:mdl.mc
            rng(seed(idx));
            [Yb,bix] = datasample(Y,nm,1,'Replace',true);% random sample with replication
            Xb = X(bix,:);
            oobix=setdiff(1:nm,bix); %% The out-of-bag indices
            Xo = X(oobix,:);Yo = Y(oobix);%% The out-of-bag data
            EC = fitcecoc(Xb,Yb,'Learners',L,'coding',coding);
            for ii=1:numel(C)
               fs(:,ii,idx) = EC.BinaryLearners{ii}.Beta;
            end
            ECNull = fitcecoc(Xb,Yb(randperm(size(Yb,1))),'Learners',L,'coding',coding);
            cvEC = crossval(EC,'kfold', k);
            %% predict
            yb = predict(EC,Xb);% assess overfitting
            yo = predict(EC,Xo);% assess oob
            yT = predict(EC,XT);% assess Test
            ycvt = cvEC.kfoldPredict;% assess crossvalidation
            yon = predict(ECNull,Xo);% assess null
            if ~isempty(XR); yR = predict(EC,XR);end% assess Replication
            %% measure perfomance 
            loss(idx) = cvEC.kfoldLoss;
            F1(idx,:) = [get.F1(Yo,yon),get.F1(Yb,yb),get.F1(Yo,yo),get.F1(Yb,ycvt),get.F1(YT,yT)];
            Acc(idx,:) = [get.acc(Yo,yon),get.acc(Yb,yb),get.acc(Yo,yo),get.acc(Yb,ycvt),get.acc(YT,yT)];
            if ~isempty(XR);Replicate(idx,:) = [get.acc(YR,yR),get.F1(YR,yR)];end           
            if ~isempty(mdl.subgroup)               
                Sb = S(bix);
                g = unique(S);
                % estimate subgroup specific models
                for ii=1:numel(g)+1
                   % for each element in subgroup train the model only by it and test on the oob & replication 
                   if ii>1
                   ix = Sb==g(ii-1);
                    ec = fitcecoc(Xb(ix,:),Yb(ix),'Learners',L,'coding',coding);
                   else
                       ec =EC;
                   end
                    for jj=1:numel(C)
                        fsg(:,jj,ii,idx) = ec.BinaryLearners{jj}.Beta;
                    end
                   for jj=1:numel(g)
                        ixt = ST==g(jj);
                        [a,b] = get.F1(YT(ixt),predict(ec,XT(ixt,:)));
                        gF1(idx,ii,jj,:) = [a,b];
                        gAcc(idx,ii,jj) = get.acc(YT(ixt),predict(ec,XT(ixt,:)));
                        if ~isempty(XR)
                           ixg = SR==g(jj);
                           [a,b] = get.F1(YR(ixg),predict(ec,XR(ixg,:)));
                           gRF1(idx,ii,jj,:) = [a,b];%mc,model,group,class 
                           gRAcc(idx,ii,jj) = get.acc(YR(ixg),predict(ec,XR(ixg,:)));
                        end
                   end
                end
            end
            if mod(idx,80);fprintf('*');else; fprintf('\n');end
        end
        toc
        output = struct('loss',loss,'F1',F1,'Acc',Acc);
        output.fs = fs;
        if ~isempty(XR);output.Replicate = Replicate;end
        if ~isempty(mdl.subgroup)
            output.gF1 = gF1;
            output.gAcc = gAcc;
            if exist('SR','var')
                output.gRF1 = gRF1;
                output.gRAcc = gRAcc;
                output.g = g;
            end
            output.fsg = fsg;
        end 
        
        save(sprintf('%s%s%s_%i.%s',partmp,filesep,mdl.name,id,'mat'),'output','-V7.3');
    otherwise
        ME = MException('Fusion:noSuchRunID','unknown run id: %s',id);
        throw(ME)
end
end


function mdl = setparams(obj)

prop = {'name','root','results','X','Y','ix','coding','Learners','kfold','ho','perm','mc','seed','stratify','YT','XT','YR','XR','subgroup'};

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
        case 'subgroup';if ia(ii);mdl.subgroup=obj.(tmp);else; mdl.subgroup=[];end
        case 'ix';if ia(ii);mdl.ix=obj.(tmp);else; mdl.ix=[];end
        case 'seed';if ia(ii);mdl.seed=obj.(tmp);else; mdl.seed=randi([1e8,1e9],mdl.mc,1);end
        case 'stratify'
            if ia(ii)&&~isempty(obj.(tmp))
                mdl.stratify=obj.(tmp);
                if mdl.perm
                    mdl.permSeed=randi([1e8,1e9],mdl.perm,1);
                    d = new.folder([],mdl.D.results,arrayfun(@(x) sprintf('perm%03d',x),1:mdl.perm,'un',0));
                    mdl.D.results = d;
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


function [X,XT,XR,Y,YT,YR,S,SR,ST,L,k,partmp,Replicate,seed,coding,nm,loss,F1,Acc,gF1,gAcc,gRF1,gRAcc]  = extractparapms(mdl)

if isfield(mdl,'permid')
    perm = mdl.permid;
    partmp = [mdl.D.results filesep sprintf('perm%03d',perm)];
    X = mdl.X(mdl.ix(:,perm),:);Y = mdl.Y(mdl.ix(:,perm));
    XT = mdl.X(~mdl.ix(:,perm),:);YT = mdl.Y(~mdl.ix(:,perm));
    if ~isempty(mdl.subgroup)
        S = mdl.subgroup.S(mdl.ix(:,perm));ST = mdl.subgroup.S(~mdl.ix(:,perm));
    else 
        S =[];ST=[];
    end
else
    partmp = [mdl.D.results, filesep,'tmp'];
    X = mdl.X(mdl.ix,:);Y = mdl.Y(mdl.ix);
    XT = mdl.X(~mdl.ix,:);YT = mdl.Y(~mdl.ix);
    if ~isempty(mdl.subgroup)
        S = mdl.subgroup.S(mdl.ix);ST = mdl.subgroup.S(~mdl.ix);
    else 
        S =[];ST=[];
    end
end

if isfield(mdl,'XR')
    XR = mdl.XR;YR = mdl.YR;
    Replicate = zeros(mdl.mc,2);
else
    XR = [];YR = [];Replicate = [];
end
C = numel(unique(Y));
L = mdl.L;k = mdl.k;
seed = mdl.seed;

[nm,~] = size(X);
coding = mdl.coding;
loss = zeros(mdl.mc,1);
F1 = zeros(mdl.mc,5);
Acc = zeros(mdl.mc,5);

if ~isempty(mdl.subgroup)
    g = unique(mdl.subgroup.S);
    gF1 = zeros(mdl.mc,numel(g),numel(g),C+1);
    gAcc = zeros(mdl.mc,numel(g),numel(g),C+1);
    if isfield(mdl.subgroup,'SR')
        gRF1 = zeros(mdl.mc,numel(g),numel(g),C+1);
        gRAcc = zeros(mdl.mc,numel(g),numel(g),C+1);
        SR = mdl.subgroup.SR;
    end
else
    gF1 = [];gAcc=[];gRF1=[];gRAcc=[];SR=[];
end


end
%------------- END OF CODE --------------
