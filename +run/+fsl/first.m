function first(T1)
%% RUN.FSL.FIRST: wrapps fsl first subcortical segmentation 
%
%    __           _             
%   / _|         (_)            
%  | |_ _   _ ___ _  ___  _ __    
%  |  _| | | / __| |/ _ \| '_ \    :- Functional and Structural 
%  | | | |_| \__ \ | (_) | | | |      Integration of Neuroimages
%  |_|  \__,_|___/_|\___/|_| |_|
%
%
%% AUTHOR:  Eyal Soreq
%  EMAIL:  e.soreq14@imperial.ac.uk
%  AFFILIATION:  Imperial College London
%  VERSION:  0.0 CREATED: 24-May-2017 13:14:23
%
%% INPUTS:
%    fn - a T1 image to run first on  
%    obj - 
%       m = 'auto','fast','none' or threshold value    default = 'auto' 
%       b = [0,1] input is already brain extracted default = 1
%       s =  a cell array relating to subcortical structures default = {'L_Accu' 'L_Amyg' 'L_Caud'  'L_Pall' 'L_Puta' 'L_Thal' 'R_Accu' 'R_Amyg' 'R_Caud'  'R_Pall' 'R_Puta' 'R_Thal','L_Hipp' 'R_Hipp'};
%
%% OUTPUT:
%
%% EXAMPLES:
%{
run.fsl.first(id,obj)
%}
%
%% DEPENDENCIES: FSL 
%
% This file is part of C^3NL Pipeline
% C^3NL Pipeline is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% C^3NL Pipeline is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with C^3NL Pipeline.If not, see <http://www.gnu.org/licenses/>.
%------------- BEGIN CODE --------------
%

[root,name] = fileparts(T1);
cd (root);
[~,fsldir] = system('echo $FSLDIR');
pth =[root,filesep,'tmp'];
c3nl.copy(T1,pth);

system(sprintf('first_flirt %s/%s %s/%s_to_std_sub',pth,name,pth,name));

roi = {'L_Accu' 'L_Amyg' 'L_Caud'  'L_Pall' 'L_Puta' 'L_Thal' 'L_Hipp' 'R_Accu' 'R_Amyg' 'R_Caud'  'R_Pall' 'R_Puta' 'R_Thal' 'R_Hipp'};
fn = c3nl.select('pth',[fsldir(1:end-1),'/data/first/models_336_bin'],'name','*.bmv');
fn = fn(~c3nl.strDetect(fn,'intref'));
n = [50 50 30 40 40 40 30 50 50 30 40 40 40 30];
for id=1:numel(roi)
scm = fn(c3nl.strDetect(fn,roi{id}));
system(sprintf('run_first -i %s/%s -t %s/%s_to_std_sub.mat -n %i -o %s/%s-%s_first -m %s',...
                            pth,name,pth,name,n(id),pth,name,roi{id},scm{1}));
end


system(sprintf('fslmerge -t %s/%s_all_fast_origsegs %s',pth,name,join(c3nl.select('name','*first*nii*'),' ')));

[Y,G] = load.vol(sprintf('%s/%s_all_fast_origsegs.nii.gz',pth,name));
y = (sum(Y,4)<100&sum(Y,4)>0).*sum(Y,4);
save.vol(y,G,sprintf('%s/%s_first.nii',root,name),'nii','uint8')


end



%------------- END OF CODE --------------
