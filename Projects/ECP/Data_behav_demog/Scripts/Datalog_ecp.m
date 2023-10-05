%% The ECP project

% MEG-ECP datalog
% Writtern by MCW group, Youssofzadeh, Vahab <vyoussofzadeh@mcw.edu>
% Lastest update: 10/5/2023

%%
clc, clear

addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/FT_fucntions/helper')

ft_path ='/opt/matlab_toolboxes/ft_packages/Stable_version/fieldtrip-master';
addpath(ft_path);
ft_defaults

indir = '/group/jbinder/ECP/MEG/MEG_Work';

cd(indir)
D = dir('./EC*');

clc
ft_progress('init', 'text',     'please wait ...');

Subinfo = [];
for i=1:length(D)
    
    ft_progress(i/length(D), 'Processing event %d from %d', i, length(D));
    pause(0.1);
    cd(fullfile(D(i).folder,D(i).name))
    Subinfo.D(i) = D(i);
    Subinfo.ID{i,:} = D(i).name;
    d_SD = rdir(['./**/tSSS','/*SD*fif']);
    d_PN = rdir(['./**/tSSS','/*PN*fif']);
    d_SM = rdir(['./**/tSSS','/*SM*fif']);
    d_anat = rdir(['./**/Anatomy','/*/brain*gz']);
    d_restEO = rdir(['./**/tSSS','/ec*RestEO*fif']);
    d_restEC = rdir(['./**/tSSS','/ec*RestEC*fif']);
    d_emproom = rdir(['./**/','/*ERNoise*fif']);
    
    if isempty(d_SD), Subinfo.task(i,1) = 0; else, Subinfo.task(i,1) = 1; end
    if isempty(d_SD), Subinfo.task(i,1) = 0; else, Subinfo.task(i,1) = 1; end
    if isempty(d_PN), Subinfo.task(i,2) = 0; else, Subinfo.task(i,2) = 1; end
    if isempty(d_SM), Subinfo.task(i,3) = 0; else, Subinfo.task(i,3) = 1; end
    if isempty(d_anat), Subinfo.anat(i,1) = 0; else, Subinfo.anat(i,1) = 1; end
    if isempty(d_restEO), Subinfo.rest(i,1) = 0; else, Subinfo.rest(i,1) = 1; end
    if isempty(d_restEC), Subinfo.rest(i,2) = 0; else, Subinfo.rest(i,2) = 1; end
    if isempty(d_emproom), Subinfo.emproom(i,1) = 0; else, Subinfo.emproom(i,1) = 1; end
    
    %- SD
    idx = 4; for ii=1:length(d_SD) if contains(lower(d_SD(ii).name), 'run1'), Subinfo.task(i,idx) = 1; break, else,Subinfo.task(i,idx) = 0; end, end
    idx = 5; for ii=1:length(d_SD) if contains(lower(d_SD(ii).name), 'run2'), Subinfo.task(i,idx) = 1; break, else,Subinfo.task(i,idx) = 0; end, end
    
    % PN
    idx = 6; for ii=1:length(d_PN), if contains(lower(d_PN(ii).name), 'run1'), Subinfo.task(i,idx) = 1; break, else,Subinfo.task(i,idx) = 0; end, end
    idx = 7; for ii=1:length(d_PN)  if contains(lower(d_PN(ii).name), 'run2'), Subinfo.task(i,idx) = 1; break, else, Subinfo.task(i,idx) = 0; end, end
    
    % SM
    idx = 8; for ii=1:length(d_SM)  if contains(lower(d_SM(ii).name), 'run1'), Subinfo.task(i,idx) = 1; break, else, Subinfo.task(i,idx) = 0; end, end
    idx = 9; for ii=1:length(d_SM)  if contains(lower(d_SM(ii).name), 'run2'), Subinfo.task(i,idx) = 1; break, else, Subinfo.task(i,idx) = 0; end, end 
    
end

%% reading subject demographic information
subdir = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/Data_behav_demog/Subj_demog/processed';

clc
[num, txt, sub_demog] = xlsread(fullfile(subdir,'ECP Progress Tracker'));

sub_list = []; k=1;
sub_cond = [];
sub_cond_idx = [];
for i=1:length(sub_demog)
    tmp = sub_demog(i,3);
    tmp2 = sub_demog(i,9);
    if ~isnan(tmp{1})
        if contains(tmp{1}, 'EC') && ~sum(isnan(tmp2{1}))
            sub_list{k} = tmp{1};
            sub_cond{k} = tmp2{1};
            if contains(sub_cond{k},'Ctrl')
                sub_cond_idx{k} = 1; sub_cond_val(k) = 1;
            elseif contains(sub_cond{k},'Patn')
                sub_cond_idx{k} = 2; sub_cond_val(k) = 2;
            end
            
            k=k+1;
        end
    end
end
disp(sub_list')
disp(sub_cond')
% disp(sub_cond_idx')

idx_ctrl = find(sub_cond_val ==1);
idx_patn = find(sub_cond_val ==2);

sub_demog_save = [sub_list',sub_cond',sub_cond_idx'];

% load('/data/MEG/Vahab/Github/MCW_MEGlab/Projects/ECP/Archive/Subj_demog/sub_demog.mat');

for i=1:length(Subinfo.ID)
    idx = find(contains(string(sub_demog_save(:,1)), Subinfo.ID{i}));
    if ~isempty(idx)
        Subinfo.group(i) = sub_demog_save(idx,2);
    end
end

%%
clc
T1 = table(Subinfo.ID); T1.Properties.VariableNames{'Var1'} = 'ID'; 
T2 = table(Subinfo.anat); T2.Properties.VariableNames{'Var1'} = 'Anat';
T3 = table(Subinfo.rest(:,1)); T3.Properties.VariableNames{'Var1'} = 'rest_EO';
T4 = table(Subinfo.rest(:,2)); T4.Properties.VariableNames{'Var1'} = 'rest_EC';
T5 = table(Subinfo.task(:,1)); T5.Properties.VariableNames{'Var1'} = 'SD';
T6 = table(Subinfo.task(:,2)); T6.Properties.VariableNames{'Var1'} = 'PN';
T7 = table(Subinfo.task(:,3)); T7.Properties.VariableNames{'Var1'} = 'SM';
T8 = table(Subinfo.emproom(:,1)); T8.Properties.VariableNames{'Var1'} = 'EmptyRoom';
T9 = table(Subinfo.group'); T9.Properties.VariableNames{'Var1'} = 'group';

% Runs 1&2
T10 = table(Subinfo.task(:,4)); T10.Properties.VariableNames{'Var1'} = 'SD_run1';
T11 = table(Subinfo.task(:,5)); T11.Properties.VariableNames{'Var1'} = 'SD_run2';
T12 = table(Subinfo.task(:,6)); T12.Properties.VariableNames{'Var1'} = 'PN_run1';
T13 = table(Subinfo.task(:,7)); T13.Properties.VariableNames{'Var1'} = 'PN_run2';
T14 = table(Subinfo.task(:,8)); T14.Properties.VariableNames{'Var1'} = 'SM_run1';
T15 = table(Subinfo.task(:,9)); T15.Properties.VariableNames{'Var1'} = 'SM_run2';


T = [T1,T9, T2,T3,T4,T5,T6,T7,T8, T10, T11, T12, T13, T14, T15];

index = [];
index.hc = find(strcmp(string(T.group), 'Ctrl'));
index.pt = find(strcmp(string(T.group), 'Patn'));
index.rest_EO = find(T.rest_EO == 1);
index.rest_EC = find(T.rest_EC == 1);
index.SD = find(T.SD == 1);
index.Anat = find(T.Anat == 1);
index.PN = find(T.PN == 1);
index.SM = find(T.SM == 1);

index.SD_1 = find(T.SD_run1 == 1);
index.SD_2 = find(T.SD_run2 == 1);

index.PN_1 = find(T.PN_run1 == 1);
index.PN_2 = find(T.PN_run2 == 1);

index.SM_1 = find(T.SM_run1 == 1);
index.SM_2 = find(T.SM_run2 == 1);

index.EmptyRoom = find(T.EmptyRoom == 1);
index.hc_anat = intersect(index.hc, index.Anat);
index.pt_anat = intersect(index.pt, index.Anat);

index.hc_anat_SD = intersect(index.hc_anat, index.SD);
index.hc_anat_SM = intersect(index.hc_anat, index.SM);
index.hc_anat_PN = intersect(index.hc_anat, index.PN);
index.hc_anat_rest_EO = intersect(index.hc_anat, index.rest_EO);
index.hc_anat_rest_EC = intersect(index.hc_anat, index.rest_EC);

index.pt_anat_SD = intersect(index.pt_anat, index.SD);
index.pt_anat_SM = intersect(index.pt_anat, index.SM);
index.pt_anat_PN = intersect(index.pt_anat, index.PN);
index.pt_anat_rest_EO = intersect(index.pt_anat, index.rest_EO);
index.pt_anat_rest_EC = intersect(index.pt_anat, index.rest_EC);
% disp(index)

% Calculate the length of each variable in the table
variable_lengths = structfun(@length, index);

% Get the field names of the structure
field_names = fieldnames(index);

% Display the variable lengths with their names
for i = 1:length(field_names)
    fprintf('%s: %d\n', field_names{i}, variable_lengths(i));
end

% disp(T.ID(index.hc_anat_SD))

%%
outdir = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/Data_behav_demog/Subj_demog';
cd(outdir)
filename = fullfile(outdir,'ECP_MEG_datalog.xlsx');
writetable(T,filename)
filename = fullfile(outdir,'ECP_MEG_datalog.csv');
disp(filename)
writetable(T,filename)



