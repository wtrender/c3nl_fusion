function [ERROR,job] = coreg(source,target,obj)
%% RUN.SPM12.COREG: runs spm12 affine co-registration between source (moving volume)
%                   to target (reference volume)  
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
%  VERSION:  0.0 CREATED: 23-May-2017 17:18:01
%
%% INPUTS:
%    source - nifti full path (if a 4D volume is selected it will treat the
%    first vol as source and the remining volumes would be automaticlly
%    defined as other (as long as other is not supplied independently
%    target - nifti full path (if a 4D volume is selected it will run coreg
%    on the first volume unless specified otherwise
%    OBJ struct  - with the following fields:
%       other - other volumes to apply the registration to 
%       targetVol - the volume to use in a 4D volume
%       sourceVol - the volume to use in a 4D volume
%       cost_fun - one of these options (default 'nmi'):
%               'nmi' - 'Normalised Mutual Information
%               'ecc' - 'Entropy Correlation Coefficeint 
%               'ncc' - 'Normalised Cross Correlation'    
%       sep - a vector representing the average sampled distance default is
%       [4 ,2] representing course sampling and than finner sampling
%       tol - The accuracy for each parameter in steps of tollarance 
%       fwhm - [s moothing over the 2D joint histogram default = [7,7]
%
%% OUTPUT:
%   ERROR - an error struct cataching spm output
%% EXAMPLES:
%{
run.spm12.coreg('source.nii','target.nii')
%}
%
%% DEPENDENCIES:
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
if ~exist('obj','var');obj=[];end
crg = getparams(source,target,obj);
jobs{1}.spm.spatial.coreg.estimate.ref = crg.ref;%{'01GI.nii,1'}
jobs{1}.spm.spatial.coreg.estimate.source = crg.source;
jobs{1}.spm.spatial.coreg.estimate.other = {crg.other};
jobs{1}.spm.spatial.coreg.estimate.eoptions = struct('cost_fun',crg.cfun,'sep',crg.sep,'tol',crg.tol,'fwhm', crg.fwhm);

s = strsplit(which(['run.spm12.' mfilename]),'+');
addpath([s{1} '3rdParty' filesep 'spm12']);
spm_jobman ('initcfg');
spm_jobman ('run', jobs);



%------------- END OF CODE --------------
end


function crg = getparams(source,target,obj)

prop = {'sourceVol','targetVol','source','target','other','cost_fun','sep','tol','fwhm','ERROR','crg'};
if ~isstruct(obj);ia = zeros(numel(prop),1); else [ia,~] = ismember(prop,fieldnames(obj));end
for ii=1:length(ia)
    tmp = prop{ii};
    switch tmp
        case 'sourceVol';if ia(ii);sourceVol=obj.(tmp);else; sourceVol=1;end
        case 'targetVol';if ia(ii);targetVol=obj.(tmp);else; targetVol=1;end
        case 'source';check.ext(source,'.nii');
        case 'target';check.ext(target,'.nii');    
        case 'other';if ia(ii);other=obj.(tmp);else; other=[];end
        case 'cost_fun';if ia(ii);cost_fun=obj.(tmp);else; cost_fun='nmi';end
        case 'sep';if ia(ii);sep=obj.(tmp);else; sep=[4,2];end
        case 'tol';if ia(ii);tol=obj.(tmp);else; tol=[0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];end
        case 'fwhm';if ia(ii);fwhm=obj.(tmp);else; fwhm=[7,7];end            
        case 'ERROR'
            % check for common problems and catch errors 
            
        case 'crg'
            % create co-reg structure 
            crg.ref = cellstr([target,',',num2str(targetVol)]);
            crg.source = cellstr([source,',',num2str(sourceVol)]);
            crg.cfun = cost_fun;
            crg.sep = sep;
            crg.tol = tol;
            crg.fwhm = fwhm;                       
            sg=load.vol(source,1);
            if numel(sg.d)>3
                crg.other = {convert.cell2mat(arrayfun(@(x) sprintf('%s,%i',source,x), 2:sg.d(4),'un',0)')};
            else
                crg.other = {''};
            end
    end
end


end