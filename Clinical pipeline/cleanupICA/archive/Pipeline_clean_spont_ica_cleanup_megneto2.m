
% ICA Cleanup Pipeline for Neuromag .fif Files
% Purpose: This script applies Independent Component Analysis (ICA) cleanup
% on neuromag .fif files typically used in MEG analysis.
% Author: MCW group, Vahab Youssof Zadeh <vyoussofzadeh@mcw.edu>
% Last Updated: 09/05/2023

%% Reset MATLAB environment
clear; clc; close('all');
warning('off');
clear; clc, close('all'); warning off

%% Initial setup
%- Input dir
restoredefaultpath

indir = '/MEG_data/epilepsy';
funcpath = '/MEG_data/MCW_pipeline/ICAcleanup_for_fif_files/func';
mnepath = '/usr/local/MATLAB_Tools/fieldtrip_2022/external/mne';

addpath(funcpath)

%- Adding path
cfg_init = [];
cfg_init.path_tools = '/usr/local/MATLAB_Tools';
allpath = do_init(cfg_init);

%%
disp('/MEG_data/epilepsy/xx/240612')
cd(indir)
[subjdir] = uigetdir;
cd(subjdir)

%%
dataConditions = {'Spont- Raw', 'Spont- (t)SSS', 'Other (sss)', 'raw only'};
for idx = 1:length(dataConditions)
    disp([num2str(idx), ': ', dataConditions{idx}]);
end
dcon = input('Select data saved format (eg 1 for Spont-Raw): ');

switch dcon
    case 1
        tag = 'spont';
        d = rdir([subjdir,['/*',tag,'*.fif']]);
    case 2
        tag = 'spont';
        d = rdir([subjdir,['/**/','sss','/*',tag,'*/*.fif']]);
    case 3
        tag = 'other';
        d = rdir(fullfile(subjdir,'/sss/*/*raw.fif'));
    case 4
        d = rdir([subjdir,'/*.fif']);
end

%%
clear subj datafolder datafile datafile1 data_disp
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
    disp('Enter data number to process:')
    datasel = input('');
else
    datasel = 1;
end
disp([subj, ' and,'])
disp(datafile1{datasel})
disp('was selected for the analysis.')

datafile = datafile1{datasel};

%% neuromag306mag layout
cfg = [];
cfg.layout = 'neuromag306mag.lay';
lay = ft_prepare_layout(cfg);
% ft_layoutplot(cfg);

%% Read fif file
cfg = []; cfg.channel = {'MEG'};
cfg.datafile  = datafile;
f_data  = ft_preprocessing(cfg);

%% apply ICA with n components
nIC = input('num of ICA (e.g., 10): ');

cfg = []; cfg.savepath = []; cfg.savefile = []; cfg.saveflag = 2; cfg.overwrite = 2;
cfg.lay = lay; cfg.n   = nIC; cfg.subj = subj; cfg.allpath = allpath; cfg.select = 1;
[cln_data, bic] = do_ica(cfg, f_data);

%%
if ~isempty(bic)
    
    % Export to fif format and Inspect it!
    addpath(mnepath)
    addpath(funcpath)
    
    tkz = tokenize(datafile,'/');
    str = strfind(datafile,'/'); savedir = datafile(1:str(end)-1);
    str_idx = strfind(tkz{end},'_raw');
    name = [tkz{end}(1:str_idx(end)), 'ic_raw.fif'];
    outfile = fullfile(savedir, name);
    
    disp(outfile)
    sok = input('name looking ok (Yes=1, No=0)?');
    
    if sok == 0
        disp('enter the new file name (e.g., xxx.fif):')
        name = input('','s');
        outfile = fullfile(savedir, name);
    end
    
    cfg = [];
    cfg.infile = datafile;
    cfg.outfile = outfile;
    cfg.cln_data = cln_data;
    do_mne_ex_read_write_raw(cfg);
    cd(savedir)
    disp('completed, data are ready to review in MEG_clinic!')
    disp(outfile);
    
    % Inspect cleaned data
    command = ['mbrowse ', outfile];
    system(command)
    
    %% Check the header file
    % dataheader = ft_read_header(datafile);
    
end
% Indicate the end of script execution
disp('Script execution completed!');
close all
