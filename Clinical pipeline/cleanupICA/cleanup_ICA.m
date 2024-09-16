%% Setup
clear; clc; close all;
warning('off');
restoredefaultpath;

indir = '/MEG_data/epilepsy';
funcpath = '/MEG_data/MCW_pipeline/Preprocess/func';
mnepath = '/usr/local/MATLAB_Tools/fieldtrip_2022/external/mne';
addpath(funcpath);

cfg_init = [];
cfg_init.path_tools = '/usr/local/MATLAB_Tools';
allpath = do_init(cfg_init);

disp('/MEG_data/epilepsy/xx/240612')
cd(indir);
[subjdir] = uigetdir;
cd(subjdir);

%% Flag analysis

flag = [];
flag.mnebrowse = 0;
flag.customfilename = 0;

%% Selection of data condition
dataConditions = {'raw', '(t) sss', 'all'};
for idx = 1:length(dataConditions)
    disp([num2str(idx), ': ', dataConditions{idx}]);
end
dcon = input('Select data condition (e.g., 1 for raw): ');

switch dcon
    case 1
        d = rdir([subjdir,'/**/*.fif']);
        d = d(~contains({d.name}, fullfile(subjdir, 'sss')));
    case 2
        d = rdir([subjdir,['/**/','sss','/**/*.fif']]);
    case 3
        d = rdir([subjdir,['/**/*.fif']]);
end

% Filter out files ending with '-eve.fif' and '-proj.fif'
d = d(~endsWith({d.name}, '-eve.fif'));
d = d(~endsWith({d.name}, '-proj.fif'));

confirmed = false;

while ~confirmed
    % Display all files for verification
    disp('Available files:')
    for i = 1:length(d)
        [pathstr, name, ext] = fileparts(d(i).name);
        disp([num2str(i), ': ', name]);
    end
    
    % Allow user to select specific files to process
    file_indices = input('Enter file numbers to process (e.g., [1 3 5]): ');
    
    % Confirm selected files before proceeding
    disp('Selected files for processing:')
    for i = file_indices
        [pathstr, name, ext] = fileparts(d(i).name);
        disp([num2str(i), ': ', name]);
    end
    confirm = input('Confirm files (Yes=1/No=0, Cancel=2): ');
    if confirm == 1
        confirmed = true;
    elseif confirm == 2
        disp('Operation cancelled.');
        return;
    else
        disp('Please reselect the files.');
    end
end

% Ask for the number of independent components
n_components = input('Enter the number of ICs to use (e.g., 10): ');

% Neuromag layout
cfg = []; cfg.layout = 'neuromag306mag.lay'; neuromaglayout  = ft_prepare_layout(cfg);

%% Process selected files
for i = file_indices
    datafile = d(i).name;
    disp(['Processing: ', datafile]);
    
    % Read fif file
    cfg = [];
    cfg.channel = {'MEG'};
    cfg.datafile = datafile;
    f_data = ft_preprocessing(cfg);
    
    % Apply ICA with specified number of components
    cfg = [];
    cfg.savepath = [];
    cfg.savefile = [];
    cfg.saveflag = 2;
    cfg.overwrite = 2;
    cfg.lay = neuromaglayout;
    cfg.n = n_components;
    cfg.subj = subjdir;
    cfg.allpath = allpath;
    cfg.select = 1;
    [cln_data, bic] = do_ica(cfg, f_data);
    
    if ~isempty(bic)
        % Export to fif format and inspect it
        addpath(mnepath);
        addpath(funcpath);
        
        tkz = tokenize(datafile, '/');
        str = strfind(datafile, '/');
        savedir = datafile(1:str(end)-1);
        str_idx = strfind(tkz{end}, '_raw');
        name = [tkz{end}(1:str_idx(end)), 'ic_raw.fif'];
        outfile = fullfile(savedir, name);
        
        disp(outfile);
        
        if flag.customfilename == 1
            sok = input('Name looking ok (Yes=1, No=0)?');
            
            if sok == 0
                disp('Enter the new file name (e.g., xxx.fif):');
                name = input('', 's');
                outfile = fullfile(savedir, name);
            end
        end
        
        cfg = [];
        cfg.infile = datafile;
        cfg.outfile = outfile;
        cfg.cln_data = cln_data;
        do_mne_ex_read_write_raw(cfg);
        cd(savedir);
        disp('Completed, data are ready to review in MEG_clinic!');
        disp(outfile);
        
        if flag.mnebrowse == 1
            % Inspect cleaned data
            command = ['mbrowse ', outfile];
            system(command);
        end
    end
end

% Indicate the end of script execution
disp('Script execution completed!');
close all;