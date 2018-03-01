function  perf=nullsvmecoc(id,data,niter,nperm,datatype,out)

%load('/Users/eyalsoreq/Dropbox/ImperialPhD/papers/WM/NN/Data/shen268/data_001.mat');
if isempty(id)
    G.sn = c3nl_select('pth',in,'type','f','name','data_*','outtype','cell');
    if nargin<5;N = length(G.sn);end
	if ~isfield(G,'jid');G.jid = '1111111';end
	cmd = {
	['#$ -hold_jid ' G.jid]
	['data=' data]
    ['niter=' niter]
	['nperm=' nperm]
	['out=' out]
	['datatype=' datatype]
	'matlab -nodisplay -nodesktop -nosplash -r "estimate.nullsvmecoc($SGE_TASK_ID,''$data'',$niter,$nperm,''$out'',''$datatype'');exit"'
	};
	obj = struct('name','myglmnet','wd','/home/eyalsoreq/code/C3NL_PhD', 'tmp',D.tmp,'cmd',{cmd},'output',D.output,'N',N,'queue','veryshort.q');
	G.jid = c3nl_qsub(obj);
else
    
    [X1,Y1,X2,Y2] = parseinput(id,data);
    perf.acc = zeros(niter,nperm+1);   
    [~,perf.seed] =system('echo $RANDOM');
    perf.seed = str2num(perf.seed);
    pref.datatype = datatype;
    rng(perf.seed); % for reproducibility 
    tic
    for ii=1:niter
        [~,mdl,~,ix,acc]=ml.fit.svmecoc(X1,Y1,'onevsall',0.2,5);% true model with random partitioning
        perf.acc(ii,1) = acc.T.f1;
        perf.acc(ii,2) = acc.t.f1;
        perf.acc(ii,3) = acc.cvt.f1;
        for jj=1:numel(mdl.BinaryLearners)
            perf.beta{ii,jj} = mdl.BinaryLearners{jj}.Beta;
        end
        [~,~,perf.acc(ii,4)]=plot.confusiongrad(Y2,predict(mdl,X2),categories(Y1));
        [~,~,~,~,acc]=ml.fit.svmecoc(X1,Y1(randperm(numel(Y1))),'onevsall',ix,5);% null model permutation with the above partitioning 
        perf.acc(ii,5) = acc.cvt.f1;
        fprintf('*');
        if ~mod(ii,80);fprintf('\n');end
    end
    toc
    if exist('out','var')
        filename = sprintf('%s%s%s.mat',out,filesep,datatype);
        save(filename,'perf');
    end
end
end


function [X1,Y1,X2,Y2] = parseinput(id,data)
    if ~isstruct(id)
        load(data);
    else
        X  = id.X;
        Y  = id.Y;
        gp = id.gp;
    end   
    gp = gp~=0;
    Y1 = Y(~gp);Y2 = Y(gp); % experiment 1 and 2 
    X1 = X(~gp,:);X2 = X(gp,:); % experiment 1 and 2 
end
