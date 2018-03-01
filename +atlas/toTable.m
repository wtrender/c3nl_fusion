function T = toTable(L,G,A)
%% ATLAS.TOTABLE: Intersects a 3D label map with litrature atlases to generate an informative table
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
%  VERSION:  0.0 CREATED: 08-Aug-2017 15:15:14
%
%% INPUTS:
%    L - a 3d label map matrix 
%    G - grid object of that map
%    A - cell of strings with the
%    following labels (default = all) :
%       1) WFU = Brodmann areas from Maldjian, J.A., Laurienti, P.J., Kraft, R.A., Burdette, J.H., 2003. An automated method for neuroanatomic and cytoarchitectonic atlas-based interrogation of fmri data sets. NeuroImage 19, 1233? 1239 (WFU Pickatlas, version X.xx).
%       2) Yeo7 & 3) Yeo17 = merged seven resting state networks from the following publications
%               A) Yeo2011_7Networks_MNI152_FreeSurferConformed1mm_LiberalMask.nii
%               B) Buckner2011_7Networks_MNI152_FreeSurferConformed1mm_LooseMask.nii
%               C) Choi2012_7Networks_MNI152_FreeSurferConformed1mm_LooseMask.nii    
%       4) AAL2 = Implementation of a new parcellation of the orbitofrontal cortex in the automated anatomical labeling atlas. Rolls ET, Joliot M & Tzourio-Mazoyer N (2015) . NeuroImage http://dx.doi.org/10.1016/j.neuroimage.2015.07.075
%       5) Shirer2012 = Shirer WR, Ryali S, Rykhlevskaia E, Menon V, Greicius MD: Decoding subject-driven cognitive states with whole-brain connectivity patterns. Cereb Cortex (2012). 
%       6) Cerb = Diedrichsen J., Maderwald S., Küper M., Thürling M., Rabe K., Gizewski ER, Ladd M, Timmann D (2011). Imaging the deep cerebellar nuclei: A probabilistic atlas and normalization procedure. Neuroimage.
%       7) HO_CS = HarvardOxford-[cort/sub]-maxprob-thr0-1mm.nii , Harvard-Oxford cortical and subcortical structural atlases
%       8) Glasser = HCP_MMP  (on the TODO list)
%       9) White matter atlas (on the TODO list)
%% OUTPUT:
%    T - table containing the following information:
%      ROIid, hemisphere, lobe, MNI(X Y Z), A(1 ... K) 
%% EXAMPLES:
%{
T = atlas.toTable(L,G,{'AAL2','WFU'});
T = atlas.toTable(L,G,{'Yeo7'})
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

source = strsplit(which(['atlas.' mfilename]),'+');
load([source{1} filesep 'atlas' filesep 'atlas.mat']);
input = fieldnames(Atlas);% use all of them
if exist('A','var');input(~ismember(input,A))=[];end
% match to source grid
for ii=1:numel(input)
   Atlas.(input{ii}).Lr = apply.interpolator(struct('V', Atlas.(input{ii}).L,'G',Atlas.(input{ii}).G),...
       struct('V','','G',G),'nearest');
end
headers = [{'ROIid', 'hemisphere', 'lobe','Volume', 'X', 'Y', 'Z'},input'];

id = unique(L(~isnan(L)));
T = table();
for ii=id(2:end)'
    tmp = table();
    tmp.ROIid = ii;
    BW = L==ii;
    ix = find(L==ii);
    tmp.Volume = numel(ix);
    [x,y,z] = ind2sub(size(L),ix);
    [I,J,K]=c3nl.deal(fix([mean(x),mean(y),mean(z)]));
    if ~BW(I,J,K) % this means that the shape is irregular and that the centroid is out of the volume
        [I,J,K]=ind2sub(size(BW),ix(knnsearch([x,y,z],[I,J,K]))); % so we find the closest 
    end
    tmp = [tmp,array2table(G.convert([J,I,K]),'VariableNames',{'X', 'Y', 'Z'})];
    
    tmp.hemisphere = categorical(find([tmp.X==0,tmp.X<0,tmp.X>0]),1:3,{'Medial','Left','Right'});
    tmp1 = Atlas.AAL2.Lr(BW);
    tmp1 = Atlas.AAL2.T(Atlas.AAL2.T.ROIid==mode(tmp1(tmp1~=0)),:);
    if ~isempty(tmp1)
        tmp.lobe = categorical(floor(tmp1.regioncode/1000),1:9,{'Nan','Frontal','Insula','Limbic','Occipital','Parietal','subcorticalGM','Temporal','Cerebelum'});
    else
        tmp.lobe = 'null';
    end
    for in=1:numel(input)
        tmp1 = Atlas.(input{in}).Lr;
        tmp1 = tmp1(ix(tmp1(ix)~=0));
        if ~isempty(tmp1)
            Lid = mode(tmp1);
            Lp = nnz(tmp1==Lid)/nnz(ix);
            Lname = Atlas.(input{in}).T.ROIname(Atlas.(input{in}).T.ROIid==Lid);
        else 
            Lid = 0;Lp = 0; Lname = {'null'};
        end
       tmp = [tmp,table(Lid,Lp,Lname,'VariableNames',cellfun(@(a) sprintf('%s_%s',a,input{in}),{'id','P','name'},'un',0))];
    end 
    T = [T;tmp];
end


%------------- END OF CODE --------------
