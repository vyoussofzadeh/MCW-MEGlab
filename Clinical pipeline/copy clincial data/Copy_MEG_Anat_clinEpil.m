clear, clear, clc, close all,

% Author, Vahab Youssof Zadeh, 2022-23
% update: 11/05/22

%% List of consented subjects (clinical MEG)

registry_dir = '/MEG_data/Research_studies/MEG_in_Epilepsy_Surgery';
cd(registry_dir)
[~, ~, raw] = xlsread('PRO00015813_Consented Patient Registry_8-17-2020.xls','Sheet1');

ft_path = '/usr/local/MATLAB_Tools/fieldtrip_20190419';
addpath(ft_path);
ft_defaults

clear sub_name
k=1;
for i=2:length(raw)
    if ~isnan(raw{i,3})
        sub_name{k} = [strtrim(raw{i,3}),'_', strtrim(raw{i,4})]; k=1+k;
    end
end
sub_name = sub_name';
disp(sub_name)

addpath(genpath('/MEG_data/LAB_MEMBERS/Vahab/Github/MCW-MEGlab/FT/Clinical pipeline/copy clincial data'))

%% Listing resting-state (tsss cleaned) raw data
indir = '/MEG_data/epilepsy';
cd(indir)
clc
d = dir('');
files = dir(indir);

clear sub_name2 data_idx
k=1;
for i=1:length(files)
    fname = files(i).name;
    if ~contains(fname,'!') && ~contains(fname,'.') ...
            && ~contains(fname,'EC') && ~contains(fname,'EMU_STDBRAIN_PROTOCOL') ....
            && ~contains(fname,'TEST_Candida') && ~contains(fname,'appletsss_kaitlintsss') ...
            && ~contains(fname,'baldwin_fossOLD') && ~contains(fname,'bednar_peggyIncomplete') ...
            && ~contains(fname,'blaszyk_sari_old') && ~contains(fname,'breier_charlotteVISIT1') ...
            && ~contains(fname,'BAK')
        
        d1 = dir(fullfile(indir,fname,'/**/brainstorm_db/data'));
        if ~isempty(d1)
            
            d_sel = [];
            for ii=1:length(d1)
                if contains(lower(d1(ii).name),fname)
                    d_sel = ii;
                    break,
                end
            end
            if ~isempty(d_sel)
                sub_name2{k} = strtrim(lower(d1(d_sel).name));
                disp(i)
                disp(sub_name2{k});
                data_idx(k) = i;
                k=1+k;
            end
        end
    end
end

%%
[sort_sub_name,~]=sort(lower(sub_name)); sort_sub_name = strtrim(sort_sub_name);
[sort_sub_name2,~]=sort(lower(sub_name2)); sort_sub_name2 = strtrim(sort_sub_name2);

[C_notfound,IA] = setdiff(sort_sub_name,sort_sub_name2);
[C,ia,ib] = intersect(sort_sub_name, sort_sub_name2, 'stable');

%% Exsining copied suject files (in temp folder)
clc
savedir = '/MEG_data/Research_studies/Epil_clinial/';
cd (savedir)
files_exist = dir(savedir);

k=1;
sub_name3 = [];
for i=1:length(files_exist)
    fname = files_exist(i).name;
    if ~contains(fname,'!') && ~contains(fname,'.')
        sub_name3{k} = fname;
        k=1+k;
    end
end
disp(sub_name3')

[sort_sub_name3,~]=sort(lower(sub_name3)); sort_sub_name3 = strtrim(sort_sub_name3);

%% Unifying lists
name_list = [];
name_list.consent = sort_sub_name;
name_list.all = sort_sub_name2';
name_list.old = sort_sub_name3';
name_list.newmatched = sort(C);
name_list.notmatched = sort(C_notfound);

%%
disp('1) copy MEG files')
disp('2) copy brainstorm db files')
disp('3) copy FS files')

disp('select task'); task = input('');
switch task
    case 1       
        Run_copyrawMEG
    case 2
        Run_copy_BrainstormDB
    case 3
        Run_copy_Freefurfer
end

%% Copy to Squiggles
cd ('/MEG_data/Research_studies/Epil_clinial/')
command = ['scp -r /MEG_data/Research_studies/Epil_clinial/jones_kellyann vyoussofzadeh@squiggles.rcc.mcw.edu:/data/MEG/Clinical/MEG_clinical_consented'];
system(command)

%%