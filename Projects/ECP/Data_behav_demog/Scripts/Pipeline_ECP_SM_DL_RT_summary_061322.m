clear; clc, close('all'); warning off,

%%
% cd '/data/MEG/Vahab/Github/MCW-MEGlab/FT';
% restoredefaultpath
% cd_org = cd;
% addpath(genpath(cd_org));
% rmpath('./Failedattemps');

cd '/data/MEG/Vahab/Github/MCW-MEGlab/FT';
restoredefaultpath
cd_org = cd;
addpath(genpath(cd_org));
% rmpath('./Failedattemps');

%- Adding path
cfg_init = [];
cfg_init.path_tools = '/data/MEG/Vahab/Github/MCW-MEGlab/tools';
[allpath, atlas] = vy_init(cfg_init);

%%
bsdir = '/group/jbinder/ECP/MEG/MEG_work_BS';
cd(bsdir)

addpath('/data/MEG/Vahab/Github/MCW-MEGlab/Projects/ECP/Run')

%%
% cd_org = '/data/MEG/Vahab/Github/MCW-MEGlab/FT/ECP';
cd_org = '/data/MEG/Vahab/Github/MCW-MEGlab/Projects/ECP';
% ECP_scriptdir = '/data/MEG/Vahab/Github/MCW-MEGlab/Projects/ECP';
ECP_scriptdir = '/data/MEG/Vahab/Github/MCW-MEGlab/Projects/ECP';

%% Performance
close all,
outdir = '/data/MEG/Research/ECP/';

%%
outdir = '/data/MEG/Research/ECP/Story_Math/';

tag = 'math';
%- Reaction time (RT)
Run_ecp_RT
%- Accuracy
Run_ecp_ACC
%- Difficulty level
Run_ecp_DL

%%
tag = 'str';
%- Reaction time (RT)
Run_ecp_RT
%- Accuracy
Run_ecp_ACC
%- Difficulty level
Run_ecp_DL

%% Summary
clc
disp(['RT_Math:', num2str(nanmean(rt.math)), '+-', num2str(nanstd(rt.math))])
disp(['RT_Str:', num2str(nanmean(rt.str)), '+-', num2str(nanstd(rt.str))])

disp('===========')
disp(['ACC_Math:', num2str(nanmean(accu_runs.math)), '+-', num2str(nanstd(accu_runs.math))])
disp(['ACC_Str:', num2str(nanmean(accu_runs.str)), '+-', num2str(nanstd(accu_runs.str))])

% disp('===========')
% disp(['trl_Math:', num2str(nanmean(accu.idx_math)), '+-', num2str(nanstd(accu.idx_math))])
% disp(['trl_Str:', num2str(nanmean(accu.idx_str)), '+-', num2str(nanstd(accu.idx_str))])

disp('===========')
disp(['DL_Math:', num2str(nanmean(DL_runs.math)), '+-', num2str(nanstd(DL_runs.math))])
disp(['DL_Math_max_min:', num2str(max(DL_val.math)), ', ', num2str(min(DL_val.math))])
disp(['DL_Str:', num2str(nanmean(DL_runs.str)), '+-', num2str(nanstd(DL_runs.str))])
disp(['DL_Str_max_min:', num2str(max(DL_val.str)), ', ', num2str(min(DL_val.str))])

%% Best subjects, DL
clc
[idx, m] = find(DL_val.math == max(DL_val.math));
max(DL_val.math)
% length(DL_val.math)
% length(sub_name)
sub_name(m)

[idx, m] = find(DL_val.str == max(DL_val.str));
max(DL_val.str)
% length(DL_val.str)
% length(sub_name)
sub_name(m)

%% Ctrl vs PT
load('/data/MEG/Vahab/Github/MCW-MEGlab/Projects/ECP/behavioural_demog/Subj_demog/sub_demog.mat');

[~, ia,ib] = intersect(sub_name1,sub_demog_save(:,1));
disp(sub_demog_save(ib,1));
disp(sub_name1(ia)');
sub_cond = sub_cond_val(ib);

idx_ctrl = find(sub_cond ==1);
idx_patn = find(sub_cond ==2);

%% acc
accu_runs.math_ctrl = accu_runs.math(idx_ctrl);
accu_runs.sub_ctrl = accu_runs.sub(idx_ctrl);

accu_runs.math_PT = accu_runs.math(idx_patn);
accu_runs.sub_PT = accu_runs.sub(idx_patn);

figure,
bar(accu_runs.math_ctrl), L = length(accu_runs.math_ctrl);
set(gca,'Xtick', 1:L,'XtickLabel',accu_runs.sub_ctrl);
set(gca,'FontSize',10,'XTickLabelRotation',90);
set(gcf, 'Position', [1000   100   1500   300]);
idx_badperf = find(accu_runs.math_ctrl < 50);
disp(sub_runs_name(idx_badperf));
set(gca,'color','none');

%%
% T = cell2table(Sub);
clear T
T = table(accu_runs.sub', accu_runs.all', accu_runs.math', accu_runs.str',sub_cond');
T.Properties.VariableNames{'Var1'} = 'ID';
T.Properties.VariableNames{'Var2'} = 'Mean_acc';
T.Properties.VariableNames{'Var3'} = 'Math_acc';
T.Properties.VariableNames{'Var4'} = 'Str_acc';
T.Properties.VariableNames{'Var5'} = 'Ctrl_PT';

savedir = '/data/MEG/Vahab/Github/MCW-MEGlab/Projects/ECP/behavioural_demog/Subj_demog/';
filename = fullfile(savedir, 'ACC_ECP_meg_runs.xlsx');
writetable(T,filename)
filename = fullfile(savedir, 'ACC_ECP_meg_runs.csv');
writetable(T,filename)
cd(savedir)

%%
[~, ia,ib] = intersect(sub_name,sub_demog_save(:,1));
disp(sub_demog_save(ib,1));
disp(sub_name(ia)');
sub_cond = sub_cond_val(ib);

idx_ctrl_run = find(sub_cond ==1);
idx_patn_run = find(sub_cond ==2);

% T = cell2table(Sub);
clear T1
T1 = table(accu.sub', accu.all',accu.math',accu.str');
T1.Properties.VariableNames{'Var1'} = 'ID';
T1.Properties.VariableNames{'Var2'} = 'Mean_acc';
T1.Properties.VariableNames{'Var3'} = 'Math_acc';
T1.Properties.VariableNames{'Var4'} = 'Str_acc';
% T1.Properties.VariableNames{'Var5'} = 'Ctrl_PT';

%%
% [~,~,c]  = unique(T.Name);
% [~,b] = unique(sort(reshape(c,size(T.Name)),2),'rows');
% T1 = T;
% T = T1(b,:);
%
% %%
% savedir = '/data/MEG/Vahab/Github/MCW-MEGlab/Projects/ECP/behavioural_demog/Subj_demog/';
% filename = fullfile(savedir, 'ACC_ECP_meg.xlsx');
% writetable(T,filename)
% filename = fullfile(savedir, 'ACC_ECP_meg.csv');
% writetable(T,filename)

savedir = '/data/MEG/Vahab/Github/MCW-MEGlab/Projects/ECP/behavioural_demog/Subj_demog/';
filename = fullfile(savedir, 'ACC_ECP_meg.xlsx');
writetable(T1,filename)
filename = fullfile(savedir, 'ACC_ECP_meg.csv');
writetable(T1,filename)
cd(savedir)

%% testing  ...
clc
SM_list = dir(fullfile(outdir,'MEG_work_ft','EC*'));

rt = [];
for i=1:length(SM_list)
    sub_name{i} = SM_list(i).name;
end
sub_name1 = unique(sub_name);

%%
load('/data/MEG/Vahab/Github/MCW-MEGlab/Projects/ECP/behavioural_demog/Subj_demog/sub_demog_ctrl.mat')
% b = accu_runs.sub_ctrl';
b = sub_name1';
a = sub_demog_ctrl_save(:,1);
[~, ia,~] = intersect(a,b);

EHQ = sub_demog_ctrl_save(ia,4);
sub_ctrl = sub_demog_ctrl_save(ia,1);


b = sub_ctrl';
a = accu_runs.sub;
[~, ia,~] = intersect(a,b);
acc_str = accu_runs.str(ia);
acc_math = accu_runs.math(ia);
acc_sub = accu_runs.sub(ia)';
acc_all = [acc_str;acc_math]';
ACC = acc_sub;
for i=1:length(acc_all)
    ACC{i,2} =  acc_all(i,1);
    ACC{i,3} =  acc_all(i,2);
end

b = sub_ctrl';
a = DL_runs.sub;
[~, ia,ib] = intersect(a,b);
DL_str = DL_runs.str(ia);
DL_math = DL_runs.math(ia);
DL_sub = DL_runs.sub(ia)';
DL_all = [DL_str;DL_math]';
DL_ctrl = DL_sub;
for i=1:length(DL_all)
    DL_ctrl{i,2} =  DL_all(i,1);
    DL_ctrl{i,3} =  DL_all(i,2);
end

b = sub_ctrl';
a = rt_runs.sub;
[~, ia,ib] = intersect(a,b);
RT_str = rt_runs.str(ia);
RT_math = rt_runs.math(ia);
RT_sub = rt_runs.sub(ia)';
RT_all = [RT_str; RT_math]';
RT_ctrl = RT_sub;
for i=1:length(RT_all)
    RT_ctrl{i,2} =  RT_all(i,1);
    RT_ctrl{i,3} =  RT_all(i,2);
end

disp([num2str(nanmean(RT_str)),' +- ', num2str(nanstd(RT_str))])
disp([num2str(nanmean(RT_math)),' +- ', num2str(nanstd(RT_math))])

%%
[num, txt, sub_demog] = xlsread('/data/MEG/Vahab/Github/MCW-MEGlab/Projects/ECP/behavioural_demog/Subj_demog/update_JB_091321/ECP_MEG_SMath_HC_EPrime_JB.xlsx');

clc
JB = [];
k=1;
for i=2:32
    JB.subID{k} = sub_demog{i,1};
    
    JB.str_DL(k) = sub_demog{i,7};
    JB.str_acc(k) = sub_demog{i,8};
    JB.str_RT(k) = sub_demog{i,11};
    
    JB.math_acc(k) = sub_demog{i,18};
    JB.math_DL(k) = sub_demog{i,15};
    JB.math_RT(k) = sub_demog{i,19};
    k=1+k;
end
% disp([num2str(mean(JB.str_acc)),' +- ', num2str(std(JB.str_acc))])
% disp([num2str(mean(JB.math_acc)),' +- ', num2str(std(JB.math_acc))])

JB1 = JB;
JB1.subID(3) = [];
JB1.str_acc(3) = [];
JB1.math_acc(3) = [];
JB1.str_DL(3) = [];
JB1.math_DL(3) = [];
JB1.str_RT(3) = [];
JB1.math_RT(3) = [];

disp([num2str(mean(JB1.str_acc)),' +- ', num2str(std(JB1.str_acc))])
disp([num2str(mean(JB1.math_acc)),' +- ', num2str(std(JB1.math_acc))])

disp([num2str(mean(JB1.str_DL)),' +- ', num2str(std(JB1.str_DL))])
disp([num2str(mean(JB1.math_DL)),' +- ', num2str(std(JB1.math_DL))])

disp([num2str(mean(JB1.str_RT)),' +- ', num2str(std(JB1.str_RT))])
disp([num2str(mean(JB1.math_RT)),' +- ', num2str(std(JB1.math_RT))])

[R,P] = corrcoef(JB.str_acc, JB.str_RT);
[R,P] = corrcoef(JB.math_acc, JB.math_RT);

%%
b = JB1.subID;
a = rt_runs.sub;
[~, ia,~] = intersect(a,b);
RT_str = rt_runs.str(ia);
RT_math = rt_runs.math(ia);
RT_sub = rt_runs.sub(ia)';
RT_all = [RT_str; RT_math]';
RT_ctrl = RT_sub;
for i=1:length(RT_all)
    RT_ctrl{i,2} =  RT_all(i,1);
    RT_ctrl{i,3} =  RT_all(i,2);
end
disp([num2str(nanmean(RT_str)),' +- ', num2str(nanstd(RT_str))])
disp([num2str(nanmean(RT_math)),' +- ', num2str(nanstd(RT_math))])

b = JB1.subID;
a = accu_runs.sub;
[~, ia,~] = intersect(a,b);
acc_str = accu_runs.str(ia);
acc_math = accu_runs.math(ia);
acc_sub = accu_runs.sub(ia)';
acc_all = [acc_str;acc_math]';
ACC = acc_sub;
for i=1:length(acc_all)
    ACC{i,2} =  acc_all(i,1);
    ACC{i,3} =  acc_all(i,2);
end
disp([num2str(nanmean(acc_str)),' +- ', num2str(nanstd(acc_str))])
disp([num2str(nanmean(acc_math)),' +- ', num2str(nanstd(acc_math))])


b = JB1.subID;
a = DL_runs.sub;
[~, ia,ib] = intersect(a,b);
DL_str = DL_runs.str(ia);
DL_math = DL_runs.math(ia);
DL_sub = DL_runs.sub(ia)';
DL_all = [DL_str;DL_math]';
DL_ctrl = DL_sub;
for i=1:length(DL_all)
    DL_ctrl{i,2} =  DL_all(i,1);
    DL_ctrl{i,3} =  DL_all(i,2);
end
disp([num2str(nanmean(DL_str)),' +- ', num2str(nanstd(DL_str))])
disp([num2str(nanmean(DL_math)),' +- ', num2str(nanstd(DL_math))])

%%
clear T
T = table(RT_sub, RT_str', RT_math', acc_str', acc_math', DL_str', DL_math');
T.Properties.VariableNames{'RT_sub'} = 'ECPID';
T.Properties.VariableNames{'Var2'} = 'RT_str';
T.Properties.VariableNames{'Var3'} = 'RT_math';
T.Properties.VariableNames{'Var4'} = 'acc_str';
T.Properties.VariableNames{'Var5'} = 'acc_math';
T.Properties.VariableNames{'Var6'} = 'DL_str';
T.Properties.VariableNames{'Var7'} = 'DL_math';

savedir = '/data/MEG/Vahab/Github/MCW_MEGlab/Projects/ECP/behavioural_demog/Subj_demog';
filename = fullfile(savedir,['behav_summ_HC','.xlsx']);
writetable(T,filename)
disp(T)