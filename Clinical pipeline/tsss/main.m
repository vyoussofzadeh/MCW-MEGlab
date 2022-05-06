
%% Post-processing Freesurfer script: workbench HCP tool
% ft_path = ('/data/MEG/Vahab/Github/fieldtrip');
% if ~ft_hastoolbox('qsub',1)
%     addpath(fullfile(ft_path,'qsub'));
% end
% addpath(genpath(fullfile('/opt/workbench/')));

%%
clear; clc, close('all'); warning off

%% Initial settings
% set(0,'DefaultFigureWindowStyle','docked')
% set(0,'DefaultFigureWindowStyle','normal')

cd('/MEG_data/Vahab/Github/MCW-MEGlab/FT');
restoredefaultpath
cd_org = cd;
addpath(genpath(cd_org));

%- Input dir
indir = '/MEG_data/epilepsy';
%- Output dir
outdir = '/MEG_data/Vahab/Processed_data';

%- Adding path
cfg_init = [];
cfg_init.path_tools = '/MEG_data/Vahab/Github/tools';
[allpath, atlas] = vy_init(cfg_init);

%%
sub_dir = '/MEG_acq/cais/copy_Vahab';
cd (sub_dir)

%%
d = rdir([sub_dir,'/*/*raw.fif']);

%%
clc
clear subj datafolder datafile datafile1
for i=1:length(d)
    [pathstr, name] = fileparts(d(i).name);
    datafolder{i} = pathstr;
    datafile{i} = d(i).name;
    Index = strfind(datafile{i}, '/');
    subj = datafile{i}(Index(3)+1:Index(4)-1);
end
datafile1 = datafile';
disp(datafile1)
disp([subj, ' data was selected for the tsss analysis'])
disp('============');
pause,

%%
cd /MEG_acq

% Strings for the command
% shell_script      = '/data/MEG/Vahab/Github/MCW-MEGlab/FT/functions/vy_anatomy_postfreesurferscript.sh';
shell_script = '/neuro/bin/util/maxfilter -gui -f $run -o /MEG_acq/logopenicppa/$1_$1/$2/tsss/$run -ctc /neuro/databases/ctc/ct_sparse.fif -cal /neuro/databases/sss/sss_cal.dat -autobad off -st 10 -corr .9';
mri_dir           = anatomy_dir;
subject_dir       = subject;

% streams_anatomy_freesurfer2.sh
command = [shell_script, ' ', mri_dir, ' ', subject_dir];

system(command);

%%


system(command)


