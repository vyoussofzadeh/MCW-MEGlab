%% The Spike Detection MEG pipline

% Spike detection MEG pipline
% Writtern by MCW group, Youssofzadeh, Vahab <vyoussofzadeh@mcw.edu>
% Lastest update: 08/12/2022

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
flag.speechanalysis = 1;     % speech analysis
flag.analysis = 1;

%% Datalog (subject details)
Datalog = [];

%% Initial settings
cd '/MEG_data/LAB_MEMBERS/Vahab/Github/MCW-MEGlab/FT';
restoredefaultpath
cd_org = cd;
addpath(genpath(cd_org));

%- Input dir
indir = '/MEG_data/epilepsy';
%- Output dir
savedir = '/MEG_data/Research_studies/Epil_annotated_data/annotated_info';


%- Adding path
cfg_init = [];
cfg_init.path_tools = '/MEG_data/LAB_MEMBERS/Vahab/Github/tools';
[allpath, atlas] = vy_init(cfg_init);


%%
[~, ~, raw] = xlsread('/MEG_data/Research_studies/MEG in Epilepsy Surgery/PRO00015813_Consented Patient Registry_8-17-2020.xls','Sheet1');

clear consented_sub
k=1;
for i=2:length(raw)
    %     disp([num2str(i), ':', num2str(length(raw))])
    if ~isnan(raw{i, 3})
        consented_sub{k} = [strtrim(raw{i, 3}),'_',strtrim(raw{i, 4})];
        k=1+k;
        %         pause,
    end
end
disp(consented_sub');

%%
cd(indir)
d = rdir([indir,['/*/','brainstorm_db','/**/*_spont*/dipoles_*.mat']]);
 
%%
% - Finding Dipole estimates
clear subj run sub_run
k=1;
for i=1:length(d)
    [pathstr, name] = fileparts(d(i).name);
    tkz = tokenize(pathstr,'/');
    tkz2 = tokenize(tkz{end},'_');
    st1 = strfind(pathstr,'/epilepsy/'); st2 = strfind(pathstr,'/brainstorm_db/');
    disp(pathstr(st1+10:st2-1))
    if ~isempty(find(contains(lower(consented_sub), lower(pathstr(st1+10:st2-1)))==1, 1))
        subj{k} = pathstr(st1+10:st2-1); %tkz{end-4};
        run{i} = tkz2{1};
        disp(subj{k})
        %     pause,
        sub_run{k,:} = [subj{k}, '_', run{i}];
        k=k+1;
    end
end
[sub_run_unq,IA,IC] = unique(sub_run);
disp(sub_run_unq);

%%
if exist(savedir, 'file') == 0, mkdir(savedir), end

clear D_annot;
for i=1:length(sub_run_unq)
    if ~exist(fullfile(savedir, [sub_run_unq{i}, '.mat']),'file')
        
        [pathstr, name] = fileparts(d(IA(i)).name);
        cd(pathstr)
        d2 = rdir('./data_Event*_trial*.mat');
        
        T_int = []; D_name = []; k=1;
        for j=1:length(d2)
            [pathstr, name] = fileparts(d2(j).name);
            clear D
            D = load(name);
            if isfield(D,'History') && size(D.History,1) >2  && ~isnan(mean(str2num(D.History{3, 3}))) && ~contains(name, '_band')
                id_l = size(D.History,1);
                if contains(D.History{id_l, 3},'baseline')
                    tkz3 = tokenize(D.History{1, 3},':');
                    D_fif = tkz3{2}(2:end-8);
                    data_ok = 1;
                    T_int(k,:) = str2num(D.History{3, 3});
%                     id_l = size(D.History,1);
                    idx = strfind(D.History{id_l, 3},'[');
                    bsl_val = str2num(D.History{id_l, 3}(idx+1:idx+4));
                    D_name{k} = name;  k=1+k;
                end
            else
                data_ok = 0;
            end
        end
        
        %%
%         command = ['mbrowse ', D_fif];
%         system(command)
        %%
        if data_ok==1 && exist(D_fif,'file')
            
            disp([num2str(i), '/', num2str(length(sub_run_unq))])
            %
            hdr = ft_read_header(D_fif);
            fs = hdr.orig.sfreq;
            first_samp = double(hdr.orig.raw.first_samp);
            if isempty(bsl_val), bsl_val = -300; end
            T_int_samplecorrected  = T_int - first_samp/fs - bsl_val/1e3;
            
            Anot = [];
            Anot.D_name = D_name;
            Anot.sname = sub_run_unq{i};
            Anot.filename = D_fif;
            Anot.T_tint = T_int_samplecorrected;
            
            %% Testing
            %         command = ['mbrowse ', D_fif]
            %         system(command)
            
            %%
            if ~exist(fullfile(savedir, [Anot.sname, '.mat']),'file')
                save(fullfile(savedir, [Anot.sname, '.mat']), 'Anot');
            end
            disp(Anot)
        end
    end
end
cd(savedir)

