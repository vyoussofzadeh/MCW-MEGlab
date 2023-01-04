
%% ECP Semantic decision task dataset, Medical College of Wisconsin

% Script: BS Process (preprocessing, source analysis)
% Project: ECP_SD
% Writtern by: Vahab Youssof Zadeh
% Update: 11/20/2022

clear; clc, close('all'); warning off,

%% Paths
restoredefaultpath
cd('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD')
addpath('./run')
Run_setpath
addpath('./data')
addpath('./run')
addpath(genpath('./functions'))

%%
% adding BS path
bs_path = '/opt/matlab_toolboxes/Brainstorm3_2021/brainstorm3'; %BS-2021
bs_path = '/opt/matlab_toolboxes/brainstorm3';

addpath(bs_path);
brainstorm
disp('choose DB from BS, then enter!');
pause

BS_dir = '/data/MEG/Research/ECP/Semantic_Decision/BS_database/';
BS_data_dir = fullfile(BS_dir,'data_all_subjects');
protocol = fullfile(BS_dir, 'data_all_subjects/protocol.mat');

%%
db_reload_database('current',1)
load(protocol);
Subj_bs = ProtocolSubjects.Subject;

L = length(Subj_bs);
k = 1;
clear subjs_bs
for i=1:length(Subj_bs)
    if ~contains(Subj_bs(i).Name, 'Group_analysis')
        datafile{k} = Subj_bs(i).FileName;
        subjs_bs{k} = Subj_bs(i).Name;
        k=1+k;
    end
end
unq_bs_subj = unique(subjs_bs);

%%
no_anat = {'EC1036'
    'EC1037'
    'EC1038'
    'EC1040'
    'EC1045'
    'EC1049'
    'EC1061'
    'EC1065'
    'EC1085'
    'EC1092'
    'EC1094'
    'EC1096'
    'EC1110'
    'EC1111'
    'EC1112'
    'EC1141'
    'EC1153'
    'EC1162'
    'EC1090'};

sub_all1 = setdiff(unq_bs_subj,no_anat);

%% delete short epochs
% cd(BS_data_dir)
% clc
% dd = rdir(fullfile('./EC*/*/results_*.mat'));
% 
% for ii=1:length(dd)
%     
%     cd(BS_data_dir)
%     [path, ~] = fileparts(dd(ii).name);
%     %     [ppath, ~] = fileparts(path);
%     if ~~exist(path,'dir')
%         cd(path)
%         
%         if ~contains(path,'@intra')
%             d = rdir('./*_average*.mat');
%             comment = []; tt=[];
%             clear tmv
%             for j=1:length(d)
%                 tmp = load(d(j).name);
%                 tt(j) = tmp.Time(end);
%                 if round(tt(j)) ~= 2
%                     rmv = 1; break,
%                 else
%                     rmv = 0;
%                 end
%             end
%             
%             if rmv ==1
%                 disp(path)
%                 pause
%                 rmdir(fullfile(BS_data_dir, path),'s')
%             end
%         end
%     end
% end

%% delete extra smoothed files
% cd(BS_data_dir)
% clc
% dd = rdir(fullfile('./Group_analysis/*/results_*.mat'));
% 
% for ii=1:length(dd)
%     disp([num2str(ii), '/', num2str(length(dd))])
%     cd(BS_data_dir)
%     %     [path, name] = fileparts(dd(ii).name);
%     %     cd(path)
%     tmp  = load(dd(ii).name);
%     %     d = rdir('./*- ssmooth*.mat');
%     
%     %     comment = [];
%     %     for j=1:length(d)
%     %         tmp = load(d(j).name);
%     % %         disp([num2str(j), ': ', tmp.Comment])
%     %         comment{j} = tmp.Comment(21);
%     %     end
%     disp(tmp.Comment)
%     if contains(tmp.Comment, '- ssmooth')
%         for j=1:length(d)
%             tmp = load(d(j).name);
%             disp([num2str(j), ': ', tmp.Comment])
%         end
%         %         del_sel = input('sel to delete:');
%         % pause
%         delete(dd(ii).name)
%     end
% end

%%
clc
close all,
subj_del = [];

cd(BS_data_dir)
dd = rdir(fullfile('./Group_analysis/ASubjects/results_*.mat'));
for jj=1:length(dd), disp([num2str(jj),':',dd(jj).name]); end

sFiles_name = [];
for jj=1:length(dd)
    disp([num2str(jj), '/' , num2str(length(dd))])
    tmp = load(dd(jj).name);
    if contains(tmp.Comment,'Symbol')
        sFiles_name{jj} = tmp.Comment(6:18);
    elseif contains(tmp.Comment,'Anim')
        sFiles_name{jj} = tmp.Comment(6:16);
    end
end

%%
Comment = []; k=1; kk=1;
need_correction = [];
for jj=1:length(sFiles_name)
    disp([num2str(jj), '/' , num2str(length(sFiles_name))])
    cd(BS_data_dir)
    tmp  = load(sFiles_name{jj});
    Comment{jj} = tmp.Comment;
    if length(tmp.Time) < 4000 && ~contains(Comment{jj}, 'mean')
        disp(length(tmp.Time))
        need_correction(k) = jj;
        k=k+1;
    end
    if contains(Comment{jj}, 'ssmooth_')
        Comment_sel{kk} = sFiles_name{jj};
        kk=kk+1;
        disp(Comment{jj})
%         pause,
    end
end


removing_sFiles  = sFiles_name(need_correction);
for i=1:length(removing_sFiles)
    removing_sFiles{i}
    [path, ~] = fileparts(removing_sFiles{i});
    if exist(path,'dir') == 7
        rmdir(fullfile(BS_data_dir, path),'s')
    end
end


%%
idx_31 = find(contains(Comment, '3_')==1); idx_32 = find(contains(Comment, '3 (')==1);
idx_21 = find(contains(Comment, '2_')==1); idx_22 = find(contains(Comment, '2 (')==1);

idx3 = [idx_31, idx_32];
idx2 = [idx_21, idx_22];

idx_ssmoth = find(contains(Comment, 'ssmooth')==1);
idx_3_smooth = intersect(idx3,idx_ssmoth);
idx_2_smooth = intersect(idx2,idx_ssmoth);

sFiles_3 = sFiles_name(idx_3_smooth);
sFiles_2 = sFiles_name(idx_2_smooth);

%%
L = length(sFiles_3); k = 1;
clear subjs_3
for i=1:length(sFiles_3)
%     disp([num2str(i), '/' , num2str(length(sFiles_3))])
    tmp = load(sFiles_3{i});
    subjs_3{k} = (tmp.Comment(9:14));
    disp([subjs_3{k}, ': ', num2str(i), '/' , num2str(length(sFiles_3))])
    k=1+k;
end
unq_bs_subj_3 = unique(subjs_3);

L = length(sFiles_2); k = 1;
clear subjs_2
for i=1:length(sFiles_2)
%     disp([num2str(i), '/' , num2str(length(sFiles_3))])
    tmp = load(sFiles_2{i});
    subjs_2{k} = (tmp.Comment(9:14));
    disp([subjs_2{k}, ': ', num2str(i), '/' , num2str(length(sFiles_2))])
    k=1+k;
end
unq_bs_subj_2 = unique(subjs_2);

%% Subject demog details
ECP_scriptdir = '/data/MEG/Vahab/Github/MCW_MEGlab/Projects/ECP';
load(fullfile(ECP_scriptdir,'behavioural_demog/Subj_demog/sub_demog.mat'));

k=1; ib_3 = [];
for j=1:length(subjs_3)
%     disp(j)
    [~, ia,ib] = intersect(subjs_3{j},sub_demog_save(:,1));
%     disp(ib)
%     pause
    if ~isempty(ib)
        ib_3(k) = ib; 
        k=k+1;
    end
end


k=1; ib_2 = [];
for j=1:length(subjs_2)
%     disp(j)
    [~, ia,ib] = intersect(subjs_2{j},sub_demog_save(:,1));
%     disp(ib)
%     pause
    if ~isempty(ib)
        ib_2(k) = ib; 
        k=k+1;
    end
end

% [~, ia,ib] = intersect(subjs_3,sub_demog_save(2,1));
% disp(sub_demog_save(ib,1)); 
sub_cond_3 = sub_cond_val(ib_3);
sub_cond_2 = sub_cond_val(ib_2);

idx_ctrl_3 = find(sub_cond_3 ==1);
idx_patn_3 = find(sub_cond_3 ==2);

idx_ctrl_2 = find(sub_cond_2 ==1);
idx_patn_2 = find(sub_cond_2 ==2);

sFiles_anim_patn = sFiles_3(idx_patn_3);
sFiles_symb_patn = sFiles_2(idx_patn_2); 

sFiles_anim_ctrl = sFiles_3(idx_ctrl_3);
sFiles_symb_ctrl = sFiles_2(idx_ctrl_2);

%% Inter-subject (group) averaging, 
stag_all = {'Anim', 'Symbol'};
%- ctrl, Anim
j = 1; conlab = 'Ctrl'; sFiles_in = sFiles_anim_ctrl;

%- Patn, Anim
j = 1; conlab = 'Patn'; sFiles_in = sFiles_anim_patn;

%- Ctrl, Symb
j = 2; conlab = 'Ctrl'; sFiles_in = sFiles_symb_ctrl;

%- Patn, Symb
j = 2; conlab = 'Patn'; sFiles_in = sFiles_symb_patn;

%%
OutputFiles = bst_process('CallProcess', 'process_average', sFiles_in, [], ...
    'avgtype',         1, ...  % Everything
    'avg_func',        1, ...  % Arithmetic average:  mean(x)
    'weighted',        0, ...
    'Comment', [stag_all{j},conlab], ...
    'scalenormalized', 0);

%% Contrasting data conditions
%- Math_vs_Stor, Ctrl,
% Process: Difference: A-B, abs
% bst_process('CallProcess', 'process_diff_ab', sFiles_anim_ctrl, sFiles_symb_ctrl, ...
%     'source_abs', 1);

%%
sFilesB = db_template('importfile');

% jj = 1; conlab = 'Ctrl'; sFiles_in = sFiles_math_ctrl;
jj = 1; conlab = 'Ctrl'; sFiles_in = sFiles_anim_ctrl;
jj = 2; conlab = 'Ctrl'; sFiles_in = sFiles_symb_ctrl;
% jj = 1; conlab = 'Patn'; sFiles_in = sFiles_math_patn;

% Process: FT t-test unequal fdr [1000ms]          H0:(A=B), H1:(A<>B)
bst_process('CallProcess', 'process_ft_sourcestatistics_VY', sFiles_in, sFilesB, ...
    'timewindow',     [1, 1], ...
    'scoutsel',       {}, ...
    'scoutfunc',      1, ...  % Mean
    'isabs',          0, ...
    'avgtime',        0, ...
    'randomizations', 5000, ...
    'statistictype',  1, ...  % Independent t-test
    'tail',           'two', ...  % Two-tailed
    'correctiontype', 4, ...  % fdr
    'minnbchan',      0, ...
    'Comment', ['FT ttest, ', stag_all{jj},conlab], ...
    'clusteralpha',   0.05);


conlab = 'Patn'; sFiles_A = sFiles_anim_patn; sFiles_B = sFiles_symb_patn;
conlab = 'Ctrl'; sFiles_A = sFiles_anim_ctrl; sFiles_B = sFiles_symb_ctrl;

% Process: FT t-test paired fdr [1000ms]          H0:(A=B), H1:(A<>B)
bst_process('CallProcess', 'process_ft_sourcestatistics_VY', sFiles_A, sFiles_B, ...
    'timewindow',     [1, 1], ...
    'scoutsel',       {}, ...
    'scoutfunc',      1, ...  % Mean
    'isabs',          0, ...
    'avgtime',        0, ...
    'randomizations', 10000, ...
    'statistictype',  2, ...  % Paired t-test
    'tail',           'two', ...  % Two-tailed
    'correctiontype', 5, ...  % max
    'minnbchan',      0, ...
    'Comment', ['FT ttest, ', num2str(length(sFiles_A)),' files, ',stag_all{1}, ' vs ', stag_all{2}, '_', conlab], ...
    'clusteralpha',   0.05);

%%
% ERD
cd('/data/MEG/Research/ECP/Semantic_Decision/BS_database/data_all_subjects/Group_analysis/@intra')
cd(fullfile(BS_dir,'data','Group_analysis/'));

loctag

anim_ctrol = load('results_average_221214_1257.mat'); 
Ianim_ctrol = anim_ctrol.ImageGridAmp;

anim_pt = load('results_average_221214_1246.mat'); 
Ianim_pt = anim_pt.ImageGridAmp; 


ttag = 'erd';

Source(Source > 0) = 0; Source = abs(Source);
[dest, name] = fileparts(OutputFiles.FileName);
tmp.Comment = [ttag, '_', tmp.Comment];
tmp.ImageGridAmp = Source;
savetag = [name,'_edit','.mat'];
save(savetag,'-struct', 'tmp'),
disp('Done, refresh BS!')

%%















%% Sentence recognition task
clear; clc, close('all'); warning off

%%
cd '/data/MEG/Vahab/Github/MCW-MEGlab/FT';
restoredefaultpath
cd_org = cd;
addpath(genpath(cd_org));
rmpath('./Failedattemps');

indir = '/group/jbinder/ECP/MEG/MEG_Work';
%- Output dir
outdir = '/group/jbinder/ECP/MEG/';

ECP_scriptdir = '/data/MEG/Vahab/Github/MCW-MEGlab/Projects/ECP';
ECP_datadir = '/data/MEG/Research/ECP/';

%- Adding path
cfg_init = [];
cfg_init.path_tools = '/data/MEG/Vahab/Github/MCW-MEGlab/tools';
[allpath, atlas] = vy_init(cfg_init);

%%
bsdir = '/group/jbinder/ECP/MEG/MEG_work_BS';
cd(bsdir)

%%
set(0,'DefaultFigureWindowStyle', 'normal')
bs_path = '/opt/matlab_toolboxes/brainstorm3';
addpath(bs_path);
brainstorm
disp('choose DB from BS, then enter!');
% pause,
% db_reload_database('current',1)

%% Read Data
disp('1:early_vs_late')
disp('2:early_vs_middle')
disp('3:middle_vs_late')
dsel = input(': ');
switch dsel
    case 1
        stag = 'early_vs_late';
        loctag = 'Results_Str_early_late';
    case 2
        stag = 'early_vs_middle';
        loctag = 'Results_Str_early_middle';
    case 3
        stag = 'middle_vs_late';
        loctag = 'Results_Str_middle_late';
end

%% listing source maps (individuals, tasks)
% stag = 'Beta';
cd(fullfile(bsdir,'data','Group_analysis',loctag));
dd = rdir('results_*.mat');

stag_all = {'str'};

sFiles_all = []; sub_all = [];
for jj = 1:length(stag_all)
    k = 1; sFiles = []; sub = [];
    for ii=1:length(dd)
        tmp = load(dd(ii).name);
        comm_data = tmp.Comment;
        tkz = tokenize(comm_data,'_');
        disp(comm_data)
%         pause,
        if contains(comm_data,stag_all{jj}(1:3)) && contains(comm_data,stag)
%             pause
            sFiles{k} = fullfile('Group_analysis', loctag, dd(ii).name); 
            sub{k} = tkz{1}(end-5:end);
            k=k+1;
        end
    end
    sFiles_all{jj} = sFiles;
    sub_all{jj} = sub;
end

%% Math and Story
% in1 = 1; in2 = 2;
% [~, ia,ib] = intersect(sub_all{in1},sub_all{in2});
% disp(sub_all{in2}(ib)'); 
% disp(sub_all{in1}(ia)');
% sFiles_math = sFiles_all{in1}(ia); sFiles_stor = sFiles_all{in2}(ib);
% sub_math = sub_all{in1}(ia);    sub_stor = sub_all{in2}(ib);

%% Subject demog details
% load('/data/MEG/Vahab/Github/MCW-MEGlab/Projects/ECP/Subj_demog/sub_demog.mat');
load(fullfile(ECP_scriptdir,'behavioural_demog/Subj_demog/sub_demog.mat'));

[~, ia,ib] = intersect(sub_all{1},sub_demog_save(:,1));
disp(sub_demog_save(ib,1)); 
% disp(sub_all(ia)');
sub_cond = sub_cond_val(ib);

idx_ctrl = find(sub_cond ==1);
idx_patn = find(sub_cond ==2);

% sFiles_math_patn = sFiles_math(idx_patn); 
sFiles_stor_patn = sFiles_all{1}(idx_patn); 

% sFiles_math_ctrl = sFiles_math(idx_ctrl); 
sFiles_stor_ctrl = sFiles_all{1}(idx_ctrl); 

%% Inter-subject (group) averaging, 

%- ctrl, math
% jj = 1; conlab = 'Ctrl'; sFiles_in = sFiles_math_ctrl;

%- Patn, math
% jj = 1; conlab = 'Patn'; sFiles_in = sFiles_math_patn;

%- Ctrl, stor
jj = 2; conlab = 'Ctrl'; sFiles_in = sFiles_stor_ctrl;

%- Patn, stor
jj = 2; conlab = 'Patn'; sFiles_in = sFiles_stor_patn;

%%
% cd(bsdir)
% Process: Average: Everything
OutputFiles = bst_process('CallProcess', 'process_average', sFiles_in, [], ...
    'avgtype',         1, ...  % Everything
    'avg_func',        1, ...  % Arithmetic average:  mean(x)
    'weighted',        0, ...
    'Comment', [stag_all{1},conlab], ...
    'scalenormalized', 0);
% ERD
cd(fullfile(bsdir,'data','Group_analysis/',loctag));

tmp = load(OutputFiles.FileName); Source = tmp.ImageGridAmp; ttag = 'erd';

Source(Source > 0) = 0; Source = abs(Source);
[dest, name] = fileparts(OutputFiles.FileName);
tmp.Comment = [ttag, '_', tmp.Comment];
tmp.ImageGridAmp = Source;
savetag = [name,'_edit','.mat'];
save(savetag,'-struct', 'tmp'),
disp('Done, refresh BS!')

%%
%- Grand avg
% for jj = 1:length(stag_all)
%     sFiles_in = sFiles_all{jj}; % Expected data condition
%     
%     % Process: Average: Everything
%     bst_process('CallProcess', 'process_average', sFiles_in, [], ...
%         'avgtype',         1, ...  % Everything
%         'avg_func',        1, ...  % Arithmetic average:  mean(x)
%         'weighted',        0, ...
%         'Comment', [stag_all{jj}], ...
%         'scalenormalized', 0);
% end

%% Contrasting data conditions
%- Math_vs_Stor, Ctrl,
% Process: Difference: A-B, abs
% bst_process('CallProcess', 'process_diff_ab', sFiles_stor_ctrl, sFiles_math_ctrl, ...
%     'source_abs', 1);

%% Contrasting data conditions
% % - Expected_vs_Unexpected
% in1 = 2; in2 = 3;
% 
% % Process: Difference: A-B, abs
% sFiles = bst_process('CallProcess', 'process_diff_ab', sFiles, sFiles2, ...
%     'source_abs', 1);

%% Independent sample t-test
sFilesB = db_template('importfile');

% jj = 1; conlab = 'Ctrl'; sFiles_in = sFiles_math_ctrl;
jj = 2; conlab = 'Ctrl'; sFiles_in = sFiles_stor_ctrl;
jj = 2; conlab = 'Patn'; sFiles_in = sFiles_stor_patn;
% jj = 1; conlab = 'Patn'; sFiles_in = sFiles_math_patn;

% Process: FT t-test unequal fdr [1000ms]          H0:(A=B), H1:(A<>B)
bst_process('CallProcess', 'process_ft_sourcestatistics_VY', sFiles_in, sFilesB, ...
    'timewindow',     [1, 1], ...
    'scoutsel',       {}, ...
    'scoutfunc',      1, ...  % Mean
    'isabs',          0, ...
    'avgtime',        0, ...
    'randomizations', 5000, ...
    'statistictype',  1, ...  % Independent t-test
    'tail',           'two', ...  % Two-tailed
    'correctiontype', 4, ...  % fdr
    'minnbchan',      0, ...
    'Comment', ['FT ttest, ', stag_all{1},conlab, '_', loctag], ...
    'clusteralpha',   0.05);

 
% for jj = 1:length(stag_all)
%     sFiles_in = sFiles_all{jj}; % Expected data condition
%     
%     % Process: FT t-test unequal fdr [1000ms]          H0:(A=B), H1:(A<>B)
%     bst_process('CallProcess', 'process_ft_sourcestatistics_VY', sFiles_in, sFilesB, ...
%         'timewindow',     [1, 1], ...
%         'scoutsel',       {}, ...
%         'scoutfunc',      1, ...  % Mean
%         'isabs',          0, ...
%         'avgtime',        0, ...
%         'randomizations', 5000, ...
%         'statistictype',  1, ...  % Independent t-test
%         'tail',           'two', ...  % Two-tailed
%         'correctiontype', 4, ...  % fdr
%         'minnbchan',      0, ...
%         'Comment', ['FT ttest, ', num2str(length(sFiles_in)),' files, cond: ',stag_all{jj}], ...
%         'clusteralpha',   0.05);
% end

%% Independent sample t-test (ROI, atlas-based)
cd(fullfile(bsdir,'data'));

sFilesB = db_template('importfile');
% atlas_lang = {'Language_tzourio_Mazoyer_02', {'Frontal_Sup_L', 'Frontal_Sup_R', 'Frontal_Sup_Orb_L', 'Frontal_Sup_Orb_R', 'Frontal_Mid_L', 'Frontal_Mid_R', 'Frontal_Mid_Orb_L', 'Frontal_Mid_Orb_R', 'Frontal_Inf_Oper_L', 'Frontal_Inf_Oper_R', 'Frontal_Inf_Tri_L', 'Frontal_Inf_Tri_R', 'Frontal_Inf_Orb_L', 'Frontal_Inf_Orb_R', 'Rolandic_Oper_L', 'Rolandic_Oper_R', 'Frontal_Sup_Medial_L', 'Frontal_Sup_Medial_R', 'Frontal_Mid_Orb_L_02', 'Frontal_Mid_Orb_R_02', 'Parietal_Sup_L', 'Parietal_Sup_R', 'Parietal_Inf_L', 'Parietal_Inf_R', 'SupraMarginal_L', 'SupraMarginal_R', 'Angular_L', 'Angular_R', 'Paracentral_Lobule_L', 'Paracentral_Lobule_R', 'Temporal_Sup_L', 'Temporal_Sup_R', 'Temporal_Pole_Sup_L', 'Temporal_Pole_Sup_R', 'Temporal_Mid_L', 'Temporal_Mid_R', 'Temporal_Pole_Mid_L', 'Temporal_Pole_Mid_R', 'Temporal_Inf_L', 'Temporal_Inf_R'}};
% atlas_lang = {'Desikan-Killiany', {'bankssts L', 'bankssts R', 'caudalanteriorcingulate L', 'caudalanteriorcingulate R', 'caudalmiddlefrontal L', 'caudalmiddlefrontal R', 'cuneus L', 'cuneus R', 'entorhinal L', 'entorhinal R', 'frontalpole L', 'frontalpole R', 'fusiform L', 'fusiform R', 'inferiorparietal L', 'inferiorparietal R', 'inferiortemporal L', 'inferiortemporal R', 'insula L', 'insula R', 'isthmuscingulate L', 'isthmuscingulate R', 'lateraloccipital L', 'lateraloccipital R', 'lateralorbitofrontal L', 'lateralorbitofrontal R', 'lingual L', 'lingual R', 'medialorbitofrontal L', 'medialorbitofrontal R', 'middletemporal L', 'middletemporal R', 'paracentral L', 'paracentral R', 'parahippocampal L', 'parahippocampal R', 'parsopercularis L', 'parsopercularis R', 'parsorbitalis L', 'parsorbitalis R', 'parstriangularis L', 'parstriangularis R', 'pericalcarine L', 'pericalcarine R', 'postcentral L', 'postcentral R', 'posteriorcingulate L', 'posteriorcingulate R', 'precentral L', 'precentral R', 'precuneus L', 'precuneus R', 'rostralanteriorcingulate L', 'rostralanteriorcingulate R', 'rostralmiddlefrontal L', 'rostralmiddlefrontal R', 'superiorfrontal L', 'superiorfrontal R', 'superiorparietal L', 'superiorparietal R', 'superiortemporal L', 'superiortemporal R', 'supramarginal L', 'supramarginal R', 'temporalpole L', 'temporalpole R', 'transversetemporal L', 'transversetemporal R'}}; 
% atlas_name = 'DK';
% 
% atlas_lang = {'Mindboggle', {'caudalanteriorcingulate L', 'caudalanteriorcingulate R', 'caudalmiddlefrontal L', 'caudalmiddlefrontal R', 'cuneus L', 'cuneus R', 'entorhinal L', 'entorhinal R', 'fusiform L', 'fusiform R', 'inferiorparietal L', 'inferiorparietal R', 'inferiortemporal L', 'inferiortemporal R', 'insula L', 'insula R', 'isthmuscingulate L', 'isthmuscingulate R', 'lateraloccipital L', 'lateraloccipital R', 'lateralorbitofrontal L', 'lateralorbitofrontal R', 'lingual L', 'lingual R', 'medialorbitofrontal L', 'medialorbitofrontal R', 'middletemporal L', 'middletemporal R', 'paracentral L', 'paracentral R', 'parahippocampal L', 'parahippocampal R', 'parsopercularis L', 'parsopercularis R', 'parsorbitalis L', 'parsorbitalis R', 'parstriangularis L', 'parstriangularis R', 'pericalcarine L', 'pericalcarine R', 'postcentral L', 'postcentral R', 'posteriorcingulate L', 'posteriorcingulate R', 'precentral L', 'precentral R', 'precuneus L', 'precuneus R', 'rostralanteriorcingulate L', 'rostralanteriorcingulate R', 'rostralmiddlefrontal L', 'rostralmiddlefrontal R', 'superiorfrontal L', 'superiorfrontal R', 'superiorparietal L', 'superiorparietal R', 'superiortemporal L', 'superiortemporal R', 'supramarginal L', 'supramarginal R', 'transversetemporal L', 'transversetemporal R'}}; 
% atlas_name = 'Mindboggle';

HO_atlas = {'HO_update_02',{'FP';'FP_02';'STGp';'STGp_02';'MTGa';'MTGa_02';'MTGp';'MTGp_02';'MTGtp';'MTGtp_02';'ITGa';'ITGa_02';'ITGp';'ITGp_02';'ITGtp';'ITGtp_02';'PoG';'PoG_02';'SPL';'SPL_02';'SmGa';'SmGa_02';'Ins';'Ins_02';'SmGp';'SmGp_02';'AG';'AG_02';'LOCs';'LOCs_02';'LOCi';'LOCi_02';'IcC';'IcC_02';'FMC';'FMC_02';'SMC';'SMC_02';'ScC';'ScC_02';'PcG';'PcG_02';'CGa';'CGa_02';'SFG';'SFG_02';'CGp';'CGp_02';'PcC';'PcC_02';'CC';'CC_02';'FOC';'FOC_02';'PhGa';'PhGa_02';'PaGp';'PaGp_02';'LG';'LG_02';'TFCa';'TFCa_02';'TFCp';'TFCp_02';'TOF';'TOF_02';'MFG';'MFG_02';'OFG';'OFG_02';'FOpC';'FOpC_02';'COpC';'COpC_02';'POpC';'POpC_02';'PP';'PP_02';'H1/H2';'H1/H2_02';'PT';'SccC';'OcP';'OcP_02';'IFGpt';'IFGpt_02';'IFGpo';'IFGpo_02';'PrG';'PrG_02';'TP';'TP_02';'STGa';'STGa_02'}};
atlas_name = 'HO';


for jj = 1:length(stag_all)
    sFiles_in = sFiles_all{jj}; % Expected data condition
      
    % Process: FT t-test unequal fdr [1000ms]          H0:(A=B), H1:(A<>B)
    bst_process('CallProcess', 'process_ft_sourcestatistics_VY', sFiles_in, sFilesB, ...
        'timewindow',     [1, 1], ...
        'scoutsel',       atlas_lang, ...
        'scoutfunc',      1, ...  % Mean
        'isabs',          0, ...
        'avgtime',        0, ...
        'randomizations', 10000, ...
        'statistictype',  1, ...  % Independent t-test
        'tail',           'two', ...  % Two-tailed
        'correctiontype', 4, ...  % fdr
        'minnbchan',      0, ...
        'Comment', ['FT ttest, ', num2str(length(sFiles_in)),' ', atlas_name, ' files, cond: ',stag_all{jj}], ...
        'clusteralpha',   0.05);   
end

%% Dependent sample t-test
% - Math_vs_Story

conlab = 'Patn'; sFiles_A = sFiles_math_patn; sFiles_B = sFiles_stor_patn;
conlab = 'Ctrl'; sFiles_A = sFiles_math_ctrl; sFiles_B = sFiles_stor_ctrl;

% Process: FT t-test paired fdr [1000ms]          H0:(A=B), H1:(A<>B)
bst_process('CallProcess', 'process_ft_sourcestatistics_VY', sFiles_A, sFiles_B, ...
    'timewindow',     [1, 1], ...
    'scoutsel',       {}, ...
    'scoutfunc',      1, ...  % Mean
    'isabs',          0, ...
    'avgtime',        0, ...
    'randomizations', 10000, ...
    'statistictype',  2, ...  % Paired t-test
    'tail',           'two', ...  % Two-tailed
    'correctiontype', 5, ...  % max
    'minnbchan',      0, ...
    'Comment', ['FT ttest, ', num2str(length(sFiles_A)),' files, ',stag_all{1}, ' vs ', stag_all{2}, '_', conlab], ...
    'clusteralpha',   0.05);
%     'correctiontype', 4, ...  % fdr


% bst_process('CallProcess', 'process_ft_sourcestatistics_VY', sFiles_A, sFiles_B, ...
%     'timewindow',     [1, 1], ...
%     'scoutsel',       {}, ...
%     'scoutfunc',      1, ...  % Mean
%     'isabs',          0, ...
%     'avgtime',        0, ...
%     'randomizations', 10000, ...
%     'statistictype',  2, ...  % Paired t-test
%     'tail',           'one+', ...  % Two-tailed
%     'correctiontype', 5, ...  % max
%     'minnbchan',      0, ...
%     'Comment', ['FT ttest, ', num2str(length(sFiles_A)),' files, ',stag_all{1}, ' vs ', stag_all{2}, conlab], ...
%     'clusteralpha',   0.05);
% %     'correctiontype', 4, ...  % fdr

%% Dependent sample t-test (ROI, atlas-based)
% - Math_vs_Story
% in1 = 1; in2 = 2;
% 
% [~, ia,ib] = intersect(sub_all{in1},sub_all{in2});
% disp(sub_all{in2}(ib)'); 
% disp(sub_all{in1}(ia)');
% sFiles_A = sFiles_all{in1}(ia); sFiles_B = sFiles_all{in2}(ib);

% atlas = {'scouts_tzourio-mazoyer_mod_052020_vy', {'Precentral_L', 'Precentral_R', 'Frontal_Sup_L', 'Frontal_Sup_R', 'Frontal_Sup_Orb_L', 'Frontal_Sup_Orb_R', 'Frontal_Mid_L', 'Frontal_Mid_R', 'Frontal_Mid_Orb_L', 'Frontal_Mid_Orb_R', 'Frontal_Inf_Oper_L', 'Frontal_Inf_Oper_R', 'Frontal_Inf_Tri_L', 'Frontal_Inf_Tri_R', 'Frontal_Inf_Orb_L', 'Frontal_Inf_Orb_R', 'Rolandic_Oper_L', 'Rolandic_Oper_R', 'Supp_Motor_Area_L', 'Supp_Motor_Area_R', 'Olfactory_L', 'Olfactory_R', 'Frontal_Sup_Medial_L', 'Frontal_Sup_Medial_R', 'Frontal_Mid_Orb_L_02', 'Frontal_Mid_Orb_R_02', 'Rectus_L', 'Rectus_R', 'Insula_L', 'Insula_R', 'Cingulum_Ant_L', 'Cingulum_Ant_R', 'Cingulum_Mid_L', 'Cingulum_Mid_R', 'Cingulum_Post_L', 'Cingulum_Post_R', 'Hippocampus_L', 'Hippocampus_R', 'ParaHippocampal_L', 'ParaHippocampal_R', 'Amygdala_L', 'Amygdala_R', 'Calcarine_L', 'Calcarine_R', 'Cuneus_L', 'Cuneus_R', 'Lingual_L', 'Lingual_R', 'Occipital_Sup_L', 'Occipital_Sup_R', 'Occipital_Mid_L', 'Occipital_Mid_R', 'Occipital_Inf_L', 'Occipital_Inf_R', 'Fusiform_L', 'Fusiform_R', 'Postcentral_L', 'Postcentral_R', 'Parietal_Sup_L', 'Parietal_Sup_R', 'Parietal_Inf_L', 'Parietal_Inf_R', 'SupraMarginal_L', 'SupraMarginal_R', 'Angular_L', 'Angular_R', 'Precuneus_L', 'Precuneus_R', 'Paracentral_Lobule_L', 'Paracentral_Lobule_R', 'Caudate_L', 'Caudate_R', 'Putamen_L', 'Putamen_R', 'Thalamus_L', 'Thalamus_R', 'Heschl_L', 'Heschl_R', 'Temporal_Sup_L', 'Temporal_Sup_R', 'Temporal_Pole_Sup_L', 'Temporal_Pole_Sup_R', 'Temporal_Mid_L', 'Temporal_Mid_R', 'Temporal_Pole_Mid_L', 'Temporal_Pole_Mid_R', 'Temporal_Inf_L', 'Temporal_Inf_R'}};
% atlas_lang = {'Language_tzourio_Mazoyer', {'Precentral_L', 'Precentral_R', 'Frontal_Sup_L', 'Frontal_Sup_R', 'Frontal_Sup_Orb_L', 'Frontal_Sup_Orb_R', 'Frontal_Mid_L', 'Frontal_Mid_R', 'Frontal_Mid_Orb_L', 'Frontal_Mid_Orb_R', 'Frontal_Inf_Oper_L', 'Frontal_Inf_Oper_R', 'Frontal_Inf_Tri_L', 'Frontal_Inf_Tri_R', 'Frontal_Inf_Orb_L', 'Frontal_Inf_Orb_R', 'Rolandic_Oper_L', 'Rolandic_Oper_R', 'Supp_Motor_Area_L', 'Supp_Motor_Area_R', 'Frontal_Sup_Medial_L', 'Frontal_Sup_Medial_R', 'Frontal_Mid_Orb_L_02', 'Frontal_Mid_Orb_R_02', 'Postcentral_L', 'Postcentral_R', 'Parietal_Sup_L', 'Parietal_Sup_R', 'Parietal_Inf_L', 'Parietal_Inf_R', 'SupraMarginal_L', 'SupraMarginal_R', 'Angular_L', 'Angular_R', 'Paracentral_Lobule_L', 'Paracentral_Lobule_R', 'Temporal_Sup_L', 'Temporal_Sup_R', 'Temporal_Pole_Sup_L', 'Temporal_Pole_Sup_R', 'Temporal_Mid_L', 'Temporal_Mid_R', 'Temporal_Pole_Mid_L', 'Temporal_Pole_Mid_R', 'Temporal_Inf_L', 'Temporal_Inf_R'}};
% atlas_lang = {'Language_tzourio_Mazoyer_02', {'Frontal_Sup_L', 'Frontal_Sup_R', 'Frontal_Sup_Orb_L', 'Frontal_Sup_Orb_R', 'Frontal_Mid_L', 'Frontal_Mid_R', 'Frontal_Mid_Orb_L', 'Frontal_Mid_Orb_R', 'Frontal_Inf_Oper_L', 'Frontal_Inf_Oper_R', 'Frontal_Inf_Tri_L', 'Frontal_Inf_Tri_R', 'Frontal_Inf_Orb_L', 'Frontal_Inf_Orb_R', 'Rolandic_Oper_L', 'Rolandic_Oper_R', 'Frontal_Sup_Medial_L', 'Frontal_Sup_Medial_R', 'Frontal_Mid_Orb_L_02', 'Frontal_Mid_Orb_R_02', 'Parietal_Sup_L', 'Parietal_Sup_R', 'Parietal_Inf_L', 'Parietal_Inf_R', 'SupraMarginal_L', 'SupraMarginal_R', 'Angular_L', 'Angular_R', 'Paracentral_Lobule_L', 'Paracentral_Lobule_R', 'Temporal_Sup_L', 'Temporal_Sup_R', 'Temporal_Pole_Sup_L', 'Temporal_Pole_Sup_R', 'Temporal_Mid_L', 'Temporal_Mid_R', 'Temporal_Pole_Mid_L', 'Temporal_Pole_Mid_R', 'Temporal_Inf_L', 'Temporal_Inf_R'}};
% atlas_lang = {'Desikan-Killiany', {'bankssts L', 'bankssts R', 'caudalanteriorcingulate L', 'caudalanteriorcingulate R', 'caudalmiddlefrontal L', 'caudalmiddlefrontal R', 'cuneus L', 'cuneus R', 'entorhinal L', 'entorhinal R', 'frontalpole L', 'frontalpole R', 'fusiform L', 'fusiform R', 'inferiorparietal L', 'inferiorparietal R', 'inferiortemporal L', 'inferiortemporal R', 'insula L', 'insula R', 'isthmuscingulate L', 'isthmuscingulate R', 'lateraloccipital L', 'lateraloccipital R', 'lateralorbitofrontal L', 'lateralorbitofrontal R', 'lingual L', 'lingual R', 'medialorbitofrontal L', 'medialorbitofrontal R', 'middletemporal L', 'middletemporal R', 'paracentral L', 'paracentral R', 'parahippocampal L', 'parahippocampal R', 'parsopercularis L', 'parsopercularis R', 'parsorbitalis L', 'parsorbitalis R', 'parstriangularis L', 'parstriangularis R', 'pericalcarine L', 'pericalcarine R', 'postcentral L', 'postcentral R', 'posteriorcingulate L', 'posteriorcingulate R', 'precentral L', 'precentral R', 'precuneus L', 'precuneus R', 'rostralanteriorcingulate L', 'rostralanteriorcingulate R', 'rostralmiddlefrontal L', 'rostralmiddlefrontal R', 'superiorfrontal L', 'superiorfrontal R', 'superiorparietal L', 'superiorparietal R', 'superiortemporal L', 'superiortemporal R', 'supramarginal L', 'supramarginal R', 'temporalpole L', 'temporalpole R', 'transversetemporal L', 'transversetemporal R'}};
% atlas_name = 'DK';
% 
% atlas_lang = {'Mindboggle', {'caudalanteriorcingulate L', 'caudalanteriorcingulate R', 'caudalmiddlefrontal L', 'caudalmiddlefrontal R', 'cuneus L', 'cuneus R', 'entorhinal L', 'entorhinal R', 'fusiform L', 'fusiform R', 'inferiorparietal L', 'inferiorparietal R', 'inferiortemporal L', 'inferiortemporal R', 'insula L', 'insula R', 'isthmuscingulate L', 'isthmuscingulate R', 'lateraloccipital L', 'lateraloccipital R', 'lateralorbitofrontal L', 'lateralorbitofrontal R', 'lingual L', 'lingual R', 'medialorbitofrontal L', 'medialorbitofrontal R', 'middletemporal L', 'middletemporal R', 'paracentral L', 'paracentral R', 'parahippocampal L', 'parahippocampal R', 'parsopercularis L', 'parsopercularis R', 'parsorbitalis L', 'parsorbitalis R', 'parstriangularis L', 'parstriangularis R', 'pericalcarine L', 'pericalcarine R', 'postcentral L', 'postcentral R', 'posteriorcingulate L', 'posteriorcingulate R', 'precentral L', 'precentral R', 'precuneus L', 'precuneus R', 'rostralanteriorcingulate L', 'rostralanteriorcingulate R', 'rostralmiddlefrontal L', 'rostralmiddlefrontal R', 'superiorfrontal L', 'superiorfrontal R', 'superiorparietal L', 'superiorparietal R', 'superiortemporal L', 'superiortemporal R', 'supramarginal L', 'supramarginal R', 'transversetemporal L', 'transversetemporal R'}}; 
% atlas_name = 'Mindboggle';

atlas_lang = {'HO_update_02',{'FP';'FP_02';'STGp';'STGp_02';'MTGa';'MTGa_02';'MTGp';'MTGp_02';'MTGtp';'MTGtp_02';'ITGa';'ITGa_02';'ITGp';'ITGp_02';'ITGtp';'ITGtp_02';'PoG';'PoG_02';'SPL';'SPL_02';'SmGa';'SmGa_02';'Ins';'Ins_02';'SmGp';'SmGp_02';'AG';'AG_02';'LOCs';'LOCs_02';'LOCi';'LOCi_02';'IcC';'IcC_02';'FMC';'FMC_02';'SMC';'SMC_02';'ScC';'ScC_02';'PcG';'PcG_02';'CGa';'CGa_02';'SFG';'SFG_02';'CGp';'CGp_02';'PcC';'PcC_02';'CC';'CC_02';'FOC';'FOC_02';'PhGa';'PhGa_02';'PaGp';'PaGp_02';'LG';'LG_02';'TFCa';'TFCa_02';'TFCp';'TFCp_02';'TOF';'TOF_02';'MFG';'MFG_02';'OFG';'OFG_02';'FOpC';'FOpC_02';'COpC';'COpC_02';'POpC';'POpC_02';'PP';'PP_02';'H1/H2';'H1/H2_02';'PT';'SccC';'OcP';'OcP_02';'IFGpt';'IFGpt_02';'IFGpo';'IFGpo_02';'PrG';'PrG_02';'TP';'TP_02';'STGa';'STGa_02'}};
atlas_name = 'HO';


sFiles_A = sFiles_math;
sFiles_B = sFiles_stor;
% Process: FT t-test paired fdr [1000ms]          H0:(A=B), H1:(A<>B)
bst_process('CallProcess', 'process_ft_sourcestatistics_VY', sFiles_A, sFiles_B, ...
    'timewindow',     [1, 1], ...
    'scoutsel',       atlas_lang, ...
    'scoutfunc',      1, ...  % Mean
    'isabs',          0, ...
    'avgtime',        0, ...
    'randomizations', 10000, ...
    'statistictype',  2, ...  % Paired t-test
    'tail',           'one-', ...  % Two-tailed
    'correctiontype', 5, ...  % max
    'minnbchan',      0, ...
    'Comment', ['FT ttest, ', num2str(length(sFiles_A)),' ', atlas_name,' files, ',stag_all{1}, ' vs ', stag_all{2}], ...
    'clusteralpha',   0.05);


conlab = 'Ctrl';
sFiles_A = sFiles_math_ctrl; 
sFiles_B = sFiles_stor_ctrl;
% Process: FT t-test paired fdr [1000ms]          H0:(A=B), H1:(A<>B)
bst_process('CallProcess', 'process_ft_sourcestatistics_VY', sFiles_A, sFiles_B, ...
    'timewindow',     [1, 1], ...
    'scoutsel',       atlas_lang, ...
    'scoutfunc',      1, ...  % Mean
    'isabs',          0, ...
    'avgtime',        0, ...
    'randomizations', 10000, ...
    'statistictype',  2, ...  % Paired t-test
    'tail',           'two', ...  % Two-tailed
    'correctiontype', 4, ...  % fdr
    'minnbchan',      0, ...
    'Comment', ['FT ttest, ', num2str(length(sFiles_A)),' ', atlas_name,' files, ',stag_all{1}, ' vs ', stag_all{2}, conlab], ...
    'clusteralpha',   0.05);

conlab = 'Patn';
sFiles_A = sFiles_math_patn; 
sFiles_B = sFiles_stor_patn;
% Process: FT t-test paired fdr [1000ms]          H0:(A=B), H1:(A<>B)
bst_process('CallProcess', 'process_ft_sourcestatistics_VY', sFiles_A, sFiles_B, ...
    'timewindow',     [1, 1], ...
    'scoutsel',       atlas_lang, ...
    'scoutfunc',      1, ...  % Mean
    'isabs',          0, ...
    'avgtime',        0, ...
    'randomizations', 10000, ...
    'statistictype',  2, ...  % Paired t-test
    'tail',           'two', ...  % Two-tailed
    'correctiontype', 5, ...  % max
    'minnbchan',      0, ...
    'Comment', ['FT ttest, ', num2str(length(sFiles_A)),' ', atlas_name,' files, ',stag_all{1}, ' vs ', stag_all{2}, conlab], ...
    'clusteralpha',   0.05);

%% plot stats on surface, color-coded map
% close all
% cd(fullfile(bsdir,'data','Group_analysis/','@intra/'));
% dd = rdir('pmatrix*.mat');
% 
% disp('1: DK, 2: Mindboggle');
% atlas_choose = input(':');
% 
% disp('1: Patn, 2: Ctrl, 3: all');
% dcon = input(':');
% 
% switch dcon
%     case 1
%         conlab = 'Patn';
%     case 2
%         conlab = 'Ctrl';
%     case 3
%         conlab = '';
% end
% 
% switch atlas_choose
%     case 1
%         % read_atlas_lang = load('/data/MEG/Vahab/Github/MCW-MEGlab/FT/attempts/Parcellation plot/scout_Language_tzourio_Mazoyer_02_42.mat');
%         read_atlas_lang = load('/data/MEG/Vahab/Github/MCW-MEGlab/FT/ECP/surface_based/scout_Desikan-Killiany_68.mat');
%         atlas_name = 'DK';
%     case 2
%         read_atlas_lang = load('/data/MEG/Vahab/Github/MCW-MEGlab/FT/ECP/surface_based/scout_Mindboggle_62.mat');
%         atlas_name = 'Mindboggle';
%     case 3
%         read_atlas_lang = load('/data/MEG/Vahab/Github/MCW-MEGlab/FT/ECP/surface_based/scout_Mindboggle_62.mat');
%         atlas_name = 'Mindboggle';
% end
% 
% nScouts = length(read_atlas_lang.Scouts);
% % colr = hsv(nScouts);
% colr = hot(nScouts);
% 
% addpath('/data/MEG/Vahab/Github/MCW-MEGlab/tools/helpful_tools/Custom colormap');
% colr = customcolormap([0 .25 .5 .75 1], {'#9d0142','#f66e45','#ffffbb','#65c0ae','#5e4f9f'},nScouts);
% colr = customcolormap(linspace(0,1,11), {'#a60126','#d7302a','#f36e43','#faac5d','#fedf8d','#fcffbf','#d7f08b','#a5d96b','#68bd60','#1a984e','#006936'},nScouts);
% colr = customcolormap(linspace(0,1,11), {'#68011d','#b5172f','#d75f4e','#f7a580','#fedbc9','#f5f9f3','#d5e2f0','#93c5dc','#4295c1','#2265ad','#062e61'},nScouts);
% clear colr
% addpath('/data/MEG/Vahab/Github/MCW-MEGlab/tools/helpful_tools');
% colr = viridis(nScouts);
% 
% thre = input('enter threshold value:');
% 
% % left and right regions
% idx_left = 1:2:nScouts;
% idx_right = 2:2:nScouts;
% idx_all = [idx_left, idx_right];
% 
% savedir = '/data/MEG/Vahab/Github/MCW-MEGlab/FT/ECP/surface_based/output_saved';
% 
% tkz = [];
% for ii=1:length(dd)
%     k=1;
%     tmp = load(dd(ii).name); L = length(tmp.tmap);
%     if L == nScouts && contains(tmp.Comment, conlab)
%         label = tmp.Description;
%         for iii=1:length(label)
%             idx1 = strfind(label{iii},'_');
%             for j=1:length(idx1)
%                 label{iii}(idx1(j)) = '-';
%             end
%         end
%         
%         disp(ii),
%         disp(tmp.Comment)
%         figure,
%         bar(abs(tmp.tmap(idx_all)))
%         set(gca,'Xtick', 1:L,'XtickLabel',label(idx_all));
%         set(gca,'FontSize',10,'XTickLabelRotation',90)
%         set(gcf, 'Position', [800   500   1500  500]);
%         grid on
%         set(gca,'color','none');
%         ylabel('t-value'); xlabel('ROI')
%         
%         [l, idx] = sort(abs(tmp.tmap),'descend');
%         read_atlas_lang1 = [];
%         for i=1:nScouts
%             if l(i) > thre*max(l)
%                 read_atlas_lang1.Scouts(k) = read_atlas_lang.Scouts(idx(i));
%                 read_atlas_lang1.Scouts(k).Color = colr(end-k+1,:);
%                 k=k+1;
% %                 disp(l(i))
%             end
%             if l(i) < thre*max(l)
%                 read_atlas_lang.Scouts(idx(i)).Color = colr(1,:);
%             end
%         end
%         
%         % -saving data for surface color-coding
%         tkz = tokenize(tmp.Comment,',');
%         read_atlas_lang1.Name = strtrim(tkz{end});
%         read_atlas_lang1.TessNbVertices = read_atlas_lang.TessNbVertices;
%         save(fullfile(savedir, ['Scout ',strtrim(tkz{end}), '_', atlas_name]),'-struct', 'read_atlas_lang1'),
%         
%         title([read_atlas_lang1.Name,conlab]);
%         saveas(gcf,fullfile(savedir, ['Scout ',strtrim(tkz{end}),'_', atlas_name,'.png']))
%     end
% end
% % cd(savedir);

%%

