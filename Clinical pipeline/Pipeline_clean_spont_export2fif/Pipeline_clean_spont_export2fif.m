clear; clc, close('all'); warning off

%% Flags
flag.preprocessing.filtering = 1;
flag.preprocessing.artifact = 1;
flag.preprocessing.ica = 1;
flag.notch = 1;
flag.freq = 1;     % TFR & FFT
flag.time = 1;     % Time-locked & Cov estimation
flag.gave = 0;     % grand average analysis
flag.anatomy = 1;     % grand average analysis
flag.sourceanalysis = 1;     % grand average analysis
flag.speechanalysis = 2;     % speech analysis

%% Initial settings
cd '/MEG_data/LAB_MEMBERS/Vahab/Github/MCW-MEGlab/FT';
restoredefaultpath
cd_org = cd;
addpath(genpath(cd_org));

%- Input dir
indir = '/MEG_data/epilepsy';

%- Adding path
cfg_init = [];
cfg_init.path_tools = '/MEG_data/LAB_MEMBERS/Vahab/Github/tools';
[allpath, atlas] = vy_init(cfg_init);

%%
cd(indir)
[subjdir] = uigetdir;
cd(subjdir)

%%
disp('1: Spont')
disp('2: SSEF')
disp('3: other')
dcon = input('sel data condition:');

switch dcon
    case 1
        tag = 'spont';
        %         d = rdir([subjdir,['/**/','sss','/*',tag,'*/*raw.fif']]);
        d = rdir([subjdir,['/**/','sss','/*',tag,'*/*.fif']]);
    case 2
        tag = 'SSEF';
        d = rdir([subjdir,['/**/','sss','/*',tag,'*/*raw*.fif']]);
    case 3
        d = rdir([subjdir,'/*sss.fif']);
end

%%
clear subj datafolder datafile datafile1
for i=1:length(d)
    [pathstr, ~] = fileparts(d(i).name);
    datafolder{i} = pathstr;
    datafile{i} = d(i).name;
    Index = strfind(datafile{i}, '/');
    subj = datafile{i}(Index(3)+1:Index(4)-1);
    data_disp{i} = [num2str(i), ': ', datafile{i}];
end
datafile1 = datafile';
disp(data_disp')
if length(datafile1) > 1
    datasel = input('choose data to analyze:');
else
    datasel = 1;
end
disp([subj, ' and,'])
disp(datafile1{datasel})
disp('was selected for the analysis.')
disp('============');

datafile = datafile1{datasel};

%%
epoch_type = 'STI101';

%% 4D layout
cfg = [];
cfg.layout = 'neuromag306mag.lay';
lay = ft_prepare_layout(cfg);
% ft_layoutplot(cfg);
disp('============');

%% ICA preprocesssing 
cfg = []; cfg.channel = {'MEG'}; 
cfg.datafile  = datafile;
% cfg.hpfreq = 0.1;
% cfg.lpfreq = 40;
f_data  = ft_preprocessing(cfg);

%%
cfg = []; cfg.savepath = []; cfg.savefile = []; cfg.saveflag = 2; cfg.overwrite = 2;
cfg.lay = lay; cfg.n   = 5; cfg.subj = subj; cfg.allpath = allpath; cfg.select = 1;
cln_data = vy_ica_cleaning_light(cfg, f_data);

%% Export to fif format
addpath('/MEG_data/LAB_MEMBERS/Vahab/Github/MCW-MEGlab/MCW_MEGlab_git/FT_fucntions/functions_new')
addpath('/usr/local/MATLAB_Tools/fieldtrip_20190419/external/mne')

tkz = tokenize(datafile,'/');
str = strfind(datafile,'/'); savedir = datafile(1:str(end)-1);
str_idx = strfind(tkz{end},'_raw');
name = [tkz{end}(1:str_idx(end)), 'ic_raw.fif'];
% name = [tkz{end}(1:end-4), '_ic.fif'];
outfile = fullfile(savedir, name);

disp(outfile)
sok = input('name looking ok (yes=1, no=0)?');

if sok == 1
    %-
    cfg = [];
    cfg.infile = datafile;
    cfg.outfile = outfile;
    cfg.cln_data = cln_data;
    do_mne_ex_read_write_raw(cfg);
    cd(savedir)
    
    disp('completed, data are ready to review in MEG_clinic!')
    disp(outfile);    
end

%% Check the header file
% dataheader = ft_read_header(datafile);
