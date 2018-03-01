function [X,Y,S,BA] = Boldbetas(roi,T,pth,zs)

%% mine.Boldbetas is a function that extracts the mean sum of betas based on the SPM flow
% Created by Eyal Soreq c3nl Imperial College London 2017
% INPUT
% =======
% roi = a 4D binary mask of different roi's
% T = a table that maps id's to nodes 
% pth = the path to a folder structure of some SPM experiment
% zs = standerdize using zscore per subject 
% OUTPUT
% =======
% X = data matrix 
% Y = response vector 
% S = subject id's
% BA = all data in a table format 


D = new.folder([],pth,{'STATS','Activity'});
if ~exist('zs','var');zs=0;end
%[L,T] = net2Dummy(Net); % load the seed roi's
d = size(roi);
n.roi = size(roi,4);
roi = reshape(roi,[],n.roi)>0;% convert roi to logical voxel x roi
rn =  matlab.lang.makeUniqueStrings(T.label);
spmFiles = c3nl.select(struct('pth',D.STATS,'name','SPM.mat'));
sn = c3nl.select(struct('pth',[pth,filesep, 'PRE'],'output','name','maxdepth',1));
BA = table();
tstart = tic;
for ii=1:numel(sn) % go over the beta images and load them
    fprintf(['-' sn{ii} '- \n' ]);
    D = new.folder(sn{ii},pth,{'STATS','Activity'});
    load(spmFiles{c3nl.strDetect(spmFiles,sn{ii})});
    xCon = struct2table(SPM.xCon);
    betaFiles = c3nl.select(struct('pth',D.STATS,'name','beta*'));
    tmp = table();
    for jj=1:height(xCon)
        ix = find(xCon.c{jj});
        vol=load.vol(betaFiles(ix,:));% load all the beta images across subjects
        info = tidyCell2Table(cellfun(@(x) strsplit(x,{' ' , '#','*'}),SPM.xX.name(ix)','un',0)); % event name
        d = size(vol);
        Y = reshape(vol,[],d(4));% reshape to voxel x samples
        y = Y;y(isnan(y))=0;
        tmp1 = array2table((y'*roi)./repmat(sum(roi),d(4),1),'variablename',matlab.lang.makeValidName(rn));
        tmp = [tmp;[info,tmp1]];
        fprintf('*');
    end
    fprintf(['finished subject ' sn{ii} '\n']);
    if zs
        temp = reshape(tmp{:,size(info,2)+1:end},[],1);
        ts = size(tmp{:,size(info,2)+1:end});
        tmp{:,size(info,2)+1:end} = reshape(zscore(temp),ts);
    end
    BA = [BA;[table(repmat(sn{ii},height(tmp),1),'variablename',{'subj'}),tmp]];
end
X = BA{:,find(c3nl.strDetect(BA.Properties.VariableNames,rn{1})):end};
Y = categorical(BA.b2);
S = categorical(cellstr(BA.subj));


end

