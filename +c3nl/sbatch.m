function jid = sbatch(obj)
%% C3NL.SBATCH: creates bash file and sumbites a job with some defined properties 
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
%  VERSION:  0.0 CREATED: 07-Jun-2017 09:06:42
%  VERSION:  0.1 updated: 09-Jun-2017 22:06:42
%
%% INPUTS:
%    obj - 
%       
%     'cmd' - the command to run as a cell array with lines
%     'workdir' - the workdir directory to store the bash file (default is tmp under pwd)
%     'jobname' - an identifier string 
%     'output' - outoput diretory to store verbos prtogress from the job (default = $PWD/output)
%     'error' - error diretory to store errors from the job (default = $PWD/output)
%     'array' - a range of job id's to run input is a string in the
%     following form : 
%               a. '0-10' (run 11 jobs starting with jid.0 and ending with jid.11
%               b. '1,4,7' (run three jobs with 1,4,7 as array id's
%               c. '1-10:2' (run five array jobs with increments of two's jid.[1,3,5,7,9])
%               d. '1-2000%50' (run 2000 array instances with a cap of 50
%               simultenous jobs at each time)
%     'time' - the time limit for the job string in the format of '01:00:00' is supplied  
%     'partition' - veryshort, short, long, verylong 
%     'depend' - job dependency a string containing previous jid
% more features to come 
%
%% OUTPUT:
%       jid = the slurm jobid 
%% EXAMPLES:
%{
obj = struct('cmd',{{'echo "My SLURM_ARRAY_TASK_ID:"$SLURM_ARRAY_TASK_ID'}},...
             'outdir','/home/eyalsoreq/slurm/',...
             'array','1,4,7','partition','veryshort')

jid = c3nl.sbatch(obj)
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
obj.origin=pwd;
slurm = setparams(obj);
[~,jid]=system(slurm.cmd);
jid = strsplit(jid,' ');
jid=jid{end};

    
end

function slurm = setparams(obj)
% Go over the different options and form a find cmd 
prop = {'dt','jobname','output','mem','cpu','error','array','time','partition','depend','cmd','workdir','newfile','obj'};
if ~isstruct(obj);ia = zeros(numel(prop),1); else; ia = ismember(prop,fieldnames(obj));end
for ii=1:length(ia)
    tmp = prop{ii};
    switch tmp
        case 'dt';dt = '#!/bin/sh';
        case 'jobname'
            if ia(ii);jn=obj.(tmp);else jn='Fusion';end 
            jobname=sprintf('#SBATCH --job-name=%s', jn);
        case 'output'
            if ia(ii);out=obj.(tmp);else;out=[pwd filesep 'output'];end
            new.folder([],out);
            output=sprintf('#SBATCH --output=%s%s%s_%%A_%%a.out', out,filesep, jn);
        case 'error'
            if ia(ii);err=obj.(tmp);else;err=out;end
            new.folder([],err);
            error=sprintf('#SBATCH --error=%s%s%s_%%A_%%a.err', err,filesep, jn);
        case 'array';if ia(ii);array=['#SBATCH --array=' obj.(tmp)];else; array='#';end    
        case 'time';if ia(ii);Time=['#SBATCH  --time=' obj.(tmp)];else; Time='#';end
        case 'partition';if ia(ii);partition=['#SBATCH --partition=' obj.(tmp)];else; partition=['#SBATCH --partition=veryshort'];end        
        case 'ntasks';ntasks=['#SBATCH --ntasks=1'];% Run on a single CPU
        case 'mem';if ia(ii);mem=['#SBATCH  --mem=' obj.(tmp)];else; mem='#SBATCH  --mem=30GB';end 
        case 'depend';if ia(ii);depend=['#SBATCH --dependency=afterok:' obj.(tmp)];else; depend='#';end        
        case 'cmd';if ia(ii);cmd=obj.(tmp);else;throw(MException('FUSION:nocommand', 'sbatch requires a command to run'));end 
        case 'cpu';if ia(ii);cpu=['#SBATCH --cpus-per-task=' num2str(obj.(tmp))];else;cpu='#SBATCH --cpus-per-task=1';end% # Number of CPU cores per task
        case 'workdir';if ia(ii);wd = ['#SBATCH  --workdir=' obj.(tmp)];else;wd = ['#SBATCH  --workdir=' pwd];end
        
        case 'tmp'
            if ia(ii);workdir = obj.(tmp);else;workdir = [pwd filesep 'tmp'];end
            new.folder([],workdir);
            cd(workdir);
        case 'newfile';[fid,fn] = new.file(jn,workdir,'w',1);
        case 'obj'            
            prop = {dt;'#';jobname;output;error;cpu;array;Time;partition;mem;wd;depend;};
            for row = 1:length(prop);fprintf(fid,'%s\n',prop{row});end
            for row = 1:length(cmd);fprintf(fid,'%s\n',cmd{row});end
            slurm.cmd = sprintf('sbatch %s',fn);
            cd(obj.origin);
            fclose(fid);
    end
end
end


%------------- END OF CODE --------------
