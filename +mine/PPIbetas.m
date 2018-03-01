function  [X,Xz,Y,S,PPI]= PPIbetas(fn,n)

%% mine.PPIbetas is a function that gets a char array of mat files loads them and spits 
% out data response vector and a table  
% Created by Eyal Soreq c3nl Imperial College London 2017
% INPUT
% =======
% fn = char matrix (each row points to a file)
% n = number of nodes in adjacency matrix 
% OUTPUT
% =======
% X = data matrix 
% Y = response vector 
% S = subject id's
% PPI = all data in a table format 


PPI = table(); % create an empty table 
ix=find(triu(true(n),1)); % construct linear indexes to the upper triangle 
[I,J]=find(triu(true(n),1)); % get their col/row index
nl = arrayfun(@(a,b) sprintf('v_%03dv_%03d',a,b),I,J,'un',0); % name them
tic
for ii=1:size(fn,1) % go over subjects
    subj = strsplit(fileparts(fn(ii,:)),'/');subj = subj{end};
    load(fn(ii,:));% load FC matrix per subj
    d = [size(Betas(1).PPI,3),numel(ix),numel(Betas)];
    cond = {};
    tmp = zeros(d(1)*d(3),d(2));% preallocate beta matrix per subj
    fprintf(['extracting subject' subj ':  ']);
    c=1;
    for jj=1:d(3) % runs
        for k=1:d(1) % trials
            temp = symAdj(Betas(jj).PPI(:,:,k),'mean');
            tmp(c,:)= temp(ix);
            c= c+1;
            fprintf('*');
        end
        fprintf('---');
        cond = [cond;Betas(jj).names];
    end
    fprintf('\n');
    subj = repmat({subj},d(3)*d(1),1);
    run = arrayfun(@(x) sprintf('S%02d',x),repelem(1:d(3),d(1))','un',0);
    cond = tidyCell2Table(cellfun(@(x) strsplit(x,{'#','_'}),cond,'un',0));
    PPI = [PPI;[table(subj),table(run),cond,array2table(tmp,'VariableNames',nl)]];
    ztmp = [ztmp;reshape(zscore(tmp(:)),size(tmp))];
end
toc
X = PPI{:,find(strDetect(PPI.Properties.VariableNames,nl{1})):end};
Xz = ztmp;
Y = categorical(PPI.b1);
S = categorical(PPI.subj);



end