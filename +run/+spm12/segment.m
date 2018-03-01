function segment(fn,obj)
%% RUN.SPM12.SEGMENT: wrapper to spm12 segmentation
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
%  VERSION:  0.0 CREATED: 24-May-2017 06:18:45
%
%% INPUTS:
%    fn - file name 
%    obj - 
%       channel - [1|2] (if you have T1 & T2 that share the same dimensions and position
%                       select 2 - this will assume that two filenames have
%                       been supplied or a 4D volume exists default = 1
%       TPM = tissue parametric maps 
%             default = SPM's TPM;  
%       biasreg - the amount of regularization applied in the segmentation
%                 process. A number needs to be selected default = .1
%       biasfwhm - the measure of bias correction to be applied default=80
%       biaswrite - save output =1 default=[0,0] (bias corrected bias
%                   field]
%       ngaus - number of gaussians to model each tissue type (the more
%               complex the more gauss needed default =[2,2,2,3,4,2]
%       mrf - the amount of cleanup to apply using markov random filed default=1
%       cleanup - the level of cleanup to be applied [0,1,2] default=1
%       reg  - steps of regularization perfomed by spm to register the tpms
%               default=[0 0.001 0.5 0.05 0.2]
%       affreg - ['','subj','mni','eastrn','none'] default='subj' which
%               means average sized brains 
%       fwhm - default=0 the recomended measure for mri
%       samp - the sampling distance smaller values are more expensive but
%              produce better results default=2
%       write - output deformation fields default=[0,0] i.e. no output  
%       native - controls the output per tissue in the native space default= [1,1,1,0,0,1]
%       dartel - controls the output per tissue in the native space default= [0,0,0,0,0,0]
%       vbm_mod    - modulated vbm in normlised space  default=[0,0,0,0,0,0]
%       vbm_umod -  modulated vbm in native space  default=[0,0,0,0,0,0]
%               
%
%% OUTPUT:
%
%% EXAMPLES:
%{
run.spm12.segment('t1_flirtCoreg.nii')
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
source = strsplit(which(['+run/+spm12/' mfilename]),'+');
addpath([source{1} '3rdParty' filesep 'spm12']);

if ~exist('obj','var');obj=[];end
seg = getparams(fn,obj);
jobs{1}.spm.spatial.preproc = struct('channel',seg.channel,'tissue',seg.tissue,'warp',seg.warp);
spm_jobman ('initcfg');
spm_jobman ('run', jobs);
%------------- END OF CODE --------------
end

function seg = getparams(fn,obj)
prop = {'channel','TPM','biasreg','biasfwhm','biaswrite','ngaus','fwhm','mrf','cleanup','reg','affreg','samp','write','native','dartel','vbm_mod','vbm_umod','seg'};
if ~isstruct(obj);ia = zeros(numel(prop),1); else [ia,~] = ismember(prop,fieldnames(obj));end
for ii=1:length(ia)
    tmp = prop{ii};
    switch tmp
        case 'channel';if ia(ii);channel=obj.(tmp);else; channel=1;end
        case 'TPM';if ia(ii);TPM=obj.(tmp);else; TPM=cellstr(spm_select('ExtFPListRec',fileparts(which('spm')),'^(TPM)*'));end
        case 'biasreg';if ia(ii);biasreg=obj.(tmp);else; biasreg=.1;end
        case 'biasfwhm';if ia(ii);biasfwhm=obj.(tmp);else; biasfwhm=80;end
        case 'biaswrite';if ia(ii);biaswrite=obj.(tmp);else; biaswrite=[0,0];end
        case 'ngaus';if ia(ii);ngaus=obj.(tmp);else; ngaus=[2,2,2,3,4,2];end
        case 'fwhm';if ia(ii);fwhm=obj.(tmp);else; fwhm=0;end                        
        case 'dartel';if ia(ii);dartel=obj.(tmp);else; dartel=[0,0,0,0,0,0];end            
        case 'mrf';if ia(ii);mrf=obj.(tmp);else; mrf=1;end            
        case 'cleanup';if ia(ii);cleanup=obj.(tmp);else; cleanup=1;end            
        case 'reg';if ia(ii);reg=obj.(tmp);else; reg=[0 0.001 0.5 0.05 0.2];end            
        case 'affreg';if ia(ii);affreg=obj.(tmp);else; affreg='subj';end            
        case 'samp';if ia(ii);samp=obj.(tmp);else; samp=2;end            
        case 'write';if ia(ii);write=obj.(tmp);else; write=[0,0];end            
        case 'native';if ia(ii);native=obj.(tmp);else; native=[1,1,1,0,0,1];end            
        case 'vbm_mod';if ia(ii);vbm_mod=obj.(tmp);else; vbm_mod=[0,0,0,0,0,0];end            
        case 'vbm_umod';if ia(ii);vbm_umod=obj.(tmp);else; vbm_umod=[0,0,0,0,0,0];end            
            
        case 'ERROR'
            % check for common problems and catch errors 
           
        case 'seg'
            % create co-reg structure 
            if channel==1
                seg.channel = struct('vols',{{fn}},'biasreg',biasreg,'biasfwhm',biasfwhm,'write',biaswrite);
            else 
                keyboard
            end
            seg.tissue(1) = struct('tpm',{TPM(1)},'ngaus',ngaus(1),'native',[native(1),dartel(1)],'warped',[vbm_mod(1),vbm_umod(1)]);
            seg.tissue(2) = struct('tpm',{TPM(2)},'ngaus',ngaus(2),'native',[native(2),dartel(2)],'warped',[vbm_mod(2),vbm_umod(2)]);
            seg.tissue(3) = struct('tpm',{TPM(3)},'ngaus',ngaus(3),'native',[native(3),dartel(3)],'warped',[vbm_mod(3),vbm_umod(3)]);
            seg.tissue(4) = struct('tpm',{TPM(4)},'ngaus',ngaus(4),'native',[native(4),dartel(4)],'warped',[vbm_mod(4),vbm_umod(4)]);
            seg.tissue(5) = struct('tpm',{TPM(5)},'ngaus',ngaus(5),'native',[native(5),dartel(5)],'warped',[vbm_mod(5),vbm_umod(5)]);
            seg.tissue(6) = struct('tpm',{TPM(6)},'ngaus',ngaus(6),'native',[native(6),dartel(6)],'warped',[vbm_mod(6),vbm_umod(6)]);
            seg.warp = struct('mrf',mrf,'cleanup',cleanup,'reg',reg,'affreg',affreg,'fwhm', fwhm,'samp',samp,'write',write);
            
    end
end


end
