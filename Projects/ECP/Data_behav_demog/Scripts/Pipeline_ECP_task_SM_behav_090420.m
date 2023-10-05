%% Story math task
clear; clc, close('all'); warning off

%%
cd '/data/MEG/Vahab/Github/MCW-MEGlab/FT';
restoredefaultpath
cd_org = cd;
addpath(genpath(cd_org));
rmpath('./Failedattemps');

indir = '/rcc/stor1/projects/ECP/MEG/MEG_Work';
%- Output dir
outdir = '/rcc/stor1/projects/ECP/MEG/';

%- Adding path
cfg_init = [];
cfg_init.path_tools = '/data/MEG/Vahab/Github/MCW-MEGlab/tools';
[allpath, atlas] = vy_init(cfg_init);

%%
bsdir = '/group/jbinder/ECP/MEG/MEG_work_BS';
cd(bsdir)

%%
% set(0,'DefaultFigureWindowStyle', 'normal')
% bs_path = '/opt/matlab_toolboxes/brainstorm3';
% addpath(bs_path);
% brainstorm
% disp('choose DB from BS, then enter!');
% pause,
% db_reload_database('current',1)

%%
% addpath('/data/MEG/Vahab/Github/MCW-MEGlab/Projects/ECP/behavioural_demog/Subj_demog')
addpath('/data/MEG/Research/ECP/Story_Math/behav')

%% listing source maps (individuals, tasks)
% cd(fullfile(bsdir,'data','Group_analysis/','@intra/'));
cd(fullfile(bsdir,'data','Group_analysis/','Beta/'));

% dd = rdir('results_average*.mat');
dd = rdir('results_*.mat');

stag_all = {'math'; 'str'};

sFiles_all = []; sub_all = [];
for jj = 1:length(stag_all)
    k = 1; sFiles = []; sub = [];
    for ii=1:length(dd)
        tmp = load(dd(ii).name);
        comm_data = tmp.Comment;
        disp(comm_data)
%         pause
        tkz = tokenize(comm_data,'/');
        if contains(comm_data,stag_all{jj}(1:3)) && contains(comm_data,'20Hz')
%             ~contains(comm_data,'files') && 
            sFiles{k} = ['Group_analysis/Beta/',dd(ii).name];
%             pause
            sub{k} = tkz{1};
            k=k+1;
        end
    end
    sFiles_all{jj} = sFiles;
    sub_all{jj} = sub;
end

%% Math and Story
in1 = 1; in2 = 2;
[~, ia,ib] = intersect(sub_all{in1},sub_all{in2});
disp(sub_all{in2}(ib)'); 
disp(sub_all{in1}(ia)');
sFiles_math = sFiles_all{in1}(ia); sFiles_stor = sFiles_all{in2}(ib);
sub_math = sub_all{in1}(ia);    sub_stor = sub_all{in2}(ib);

%% Subject demographic details
clear sub_cond_val ia ib
load('/data/MEG/Research/ECP/Story_Math/behav/sub_demog.mat');

[~, ia,ib] = intersect(sub_math,sub_demog_save(:,1));
disp(sub_demog_save(ib,1)); 
disp(sub_math(ia)');
sub_cond = sub_cond_val(ib);

idx_ctrl = find(sub_cond ==1);
idx_patn = find(sub_cond ==2);

sub_ID = sub_math(ia);

sFiles_math_patn = sFiles_math(idx_patn); 
sFiles_stor_patn = sFiles_stor(idx_patn); 

sFiles_math_ctrl = sFiles_math(idx_ctrl); 
sFiles_stor_ctrl = sFiles_stor(idx_ctrl);

sub_math_ctrl = sub_ID(idx_ctrl);
sub_math_pnts = sub_ID(idx_patn);

% disp(sub_math_ctrl')
% disp(sub_math_pnts')

%% subject demographic, control, age
clear sub_demog_ctrl_age_corrected sub_demog_ctrl_age sub_demog_ctrl_ID_corrected

load('sub_demog_ctrl.mat');
[~, ia,ib] = intersect(sub_math_ctrl,sub_demog_ctrl_save(:,1));
disp(sub_demog_ctrl_save(ib,1)); % ID
disp(sub_demog_ctrl_save(ib,3)); % age
disp(sub_math_ctrl(ia)'); % double-checking ...

sub_demog_ctrl_age = sub_demog_ctrl_save(ib,3);
sub_demog_ctrl_ID = sub_demog_ctrl_save(ib,1);

clear sub_demog_ctrl_age_corrected 
missing_data = [];
k=1;
for i=1:length(sub_demog_ctrl_age)
    if isnumeric(sub_demog_ctrl_age{i})
         sub_demog_ctrl_age_corrected(k) = sub_demog_ctrl_age{i};
         sub_demog_ctrl_ID_corrected{k}  = sub_demog_ctrl_ID{i};
         k=k+1;
    else
        missing_data = [missing_data, i]; 
    end   
end

mean(sub_demog_ctrl_age_corrected)
std(sub_demog_ctrl_age_corrected)

%% subject demographic, PT, age
clear sub_demog_PT_age_corrected sub_demog_PT_age sub_demog_PT_ID_corrected

load('sub_demog_pnts.mat');
[~, ia,ib] = intersect(sub_math_pnts,sub_demog_PT_save(:,1));
disp(sub_demog_PT_save(ib,1)); % ID
disp(sub_demog_PT_save(ib,4)); % age
disp(sub_math_pnts(ia)'); % double-checking ...

sub_demog_PT_age = sub_demog_PT_save(ib,4);
sub_demog_PT_ID = sub_demog_PT_save(ib,1);

missing_data  =[];
k=1;
for i=1:length(sub_demog_PT_age)
    if isnumeric(sub_demog_PT_age{i})
        sub_demog_PT_age_corrected(k) = sub_demog_PT_age{i};
        sub_demog_PT_ID_corrected{k}  = sub_demog_PT_ID{i};
        k=k+1;
    else
        missing_data = [missing_data, i];
    end
end

mean(sub_demog_PT_age_corrected)
std(sub_demog_PT_age_corrected)

%% Age
% cd('/data/MEG/Projects/ECP/behav')
cd('/data/MEG/Vahab/Github/MCW-MEGlab/Projects/ECP/behavioural_demog')
figure,
bar([mean(sub_demog_ctrl_age_corrected); mean(sub_demog_PT_age_corrected)]', 0.4), title('Age, Ctrl Pant')
hold on
errorbar(1:2,[mean(sub_demog_ctrl_age_corrected); mean(sub_demog_PT_age_corrected)]',[std(sub_demog_ctrl_age_corrected); std(sub_demog_PT_age_corrected)]')
set(gca,'Xtick', 1:2,'XtickLabel',{'Ctrl','Pant'});
xlabel('Sub'), ylabel('Acc (%)'),
% set(gca,'FontSize',10,'XTickLabelRotation',90);
grid on
set(gca,'color','none');
% legend({'Str';'Math'})

%% Accuracy
clc
outdir = '/data/MEG/Research/ECP/Story_Math/';
outd_sub = rdir(fullfile(outdir,'MEG_work_ft','**', 'Acc*.mat'));

accu = [];
for i=1:length(outd_sub)
    load(outd_sub(i).name);
    tkz = tokenize(outd_sub(i).name,'/');
    accu.all(i) = accuracy; accu.str(i) = accuracy_str; accu.math(i) = accuracy_math;
    sub_name{i} = tkz{end-3};
%     sub_run{i} = tkz{end-1};
end

sub_name1 = unique(sub_name);

% ave of both runs.
accu_runs = [];
for i=1:length(sub_name1)
    idx = find(strcmp(sub_name1{i},sub_name)==1);
    if length(idx) ==2
    accu_runs.all(i) = (accu.all(idx(1)) + accu.all(idx(2)))/2; 
    accu_runs.str(i) = (accu.str(idx(1)) + accu.str(idx(2)))/2; 
    accu_runs.math(i) = (accu.math(idx(1)) + accu.math(idx(2)))/2;
    tmp = sub_name(idx(1));
    sub_runs_name{i} = tmp{1};
    elseif length(idx) ==1
    accu_runs.all(i) = accu.all(idx); 
    accu_runs.str(i) = accu.str(idx); 
    accu_runs.math(i) = accu.math(idx);
    tmp = sub_name(idx);
    sub_runs_name{i} = tmp{1};
    disp(sub_name(idx))
    end
end

% Second run only.
% accu_runs = [];
% for i=1:length(sub_name1)
%     idx = find(strcmp(sub_name1{i},sub_name)==1);
%     if length(idx) ==2
%     accu_runs.all(i) = accu.all(idx(2)); 
%     accu_runs.str(i) = accu.str(idx(2)); 
%     accu_runs.math(i) = accu.math(idx(2));
%     tmp = sub_name(idx(1));
%     sub_runs_name{i} = tmp{1};
%     elseif length(idx) ==1
%     accu_runs.all(i) = accu.all(idx); 
%     accu_runs.str(i) = accu.str(idx); 
%     accu_runs.math(i) = accu.math(idx);
%     tmp = sub_name(idx);
%     sub_runs_name{i} = tmp{1};
%     disp(sub_name(idx))
%     end
% end

figure,
bar(accu_runs.all), L = length(accu_runs.all);
set(gca,'Xtick', 1:L,'XtickLabel',sub_runs_name);
set(gca,'FontSize',10,'XTickLabelRotation',90);
set(gcf, 'Position', [1000   100   1500   300]);
idx_badperf = find(accu_runs.all < 50);
disp(sub_runs_name(idx_badperf));
set(gca,'color','none');
% switch tag
%     case 'Math'
%         stag = 'math';
%     case 'Str'
%         stag = 'str';
% end
% title([tag, ', acc: ', num2str(mean(accu.(stag))), '%'])

%% Difficulty level
clc
outdir = '/data/MEG/Research/ECP/Story_Math/';
outd_sub = rdir(fullfile(outdir,'MEG_work_ft','**', 'DL*.mat'));

clear sub_name
DL_val = [];
for i=1:length(outd_sub)
    load(outd_sub(i).name);
    tkz = tokenize(outd_sub(i).name,'/');
    
%     DL_val.all(i) = DL;
    idx = find(DL_str_uniq > 4);
    if ~isempty(idx)
        disp(i)
        disp(tkz{end-3})
        DL_str_uniq(idx)=[];
    end
    DL_val.str(i) = mean(DL_str_uniq);
    DL_val.math(i) = mean(DL_math_uniq);
    
    sub_name{i} = tkz{end-3};
    %     sub_run{i} = tkz{end-1};
end

sub_name1 = unique(sub_name);

DL_runs = [];
for i=1:length(sub_name1)
    idx = find(strcmp(sub_name1{i},sub_name)==1);
    if length(idx) ==2
%     DL_runs.all(i) = (DL_val.all(idx(1)) + DL_val.all(idx(2)))/2; 
    DL_runs.str(i) = (DL_val.str(idx(1)) + DL_val.str(idx(2)))/2; 
    DL_runs.math(i) = (DL_val.math(idx(1)) + DL_val.math(idx(2)))/2;
    tmp = sub_name(idx(1));
    sub_runs_name{i} = tmp{1};
    elseif length(idx) ==1
%     DL_runs.all(i) = DL_val.all(idx); 
    DL_runs.str(i) = DL_val.str(idx); 
    DL_runs.math(i) = DL_val.math(idx);
    tmp = sub_name(idx);
    sub_runs_name{i} = tmp{1};
    disp(sub_name(idx))
    end
end

% figure,
% bar(DL_runs.str), L = length(DL_runs.str);
% set(gca,'Xtick', 1:L,'XtickLabel',sub_runs_name);
% set(gca,'FontSize',10,'XTickLabelRotation',90);
% set(gcf, 'Position', [1000   100   1500   300]);
% idx_badperf = find(DL_runs.str < 2);
% disp(sub_runs_name(idx_badperf));
% set(gca,'color','none');
% 
% figure,
% bar(DL_runs.math), L = length(DL_runs.math);
% set(gca,'Xtick', 1:L,'XtickLabel',sub_runs_name);
% set(gca,'FontSize',10,'XTickLabelRotation',90);
% set(gcf, 'Position', [1000   100   1500   300]);
% idx_badperf = find(DL_runs.math < 5);
% disp(sub_runs_name(idx_badperf));
% set(gca,'color','none');

figure,
bar((DL_runs.math + DL_runs.str)/2), L = length(DL_runs.math);
set(gca,'Xtick', 1:L,'XtickLabel',sub_runs_name);
set(gca,'FontSize',10,'XTickLabelRotation',90);
set(gcf, 'Position', [1000   100   1500   300]);
set(gca,'color','none');

%%
[~, ia,ib] = intersect(sub_ID,sub_runs_name);

sub_ID1 = sub_ID(ia);
Acc_str = accu_runs.str(ib);
Acc_math = accu_runs.math(ib);

Acc_math_patn = Acc_math(idx_patn); 
Acc_str_patn = Acc_str(idx_patn); 

Acc_math_ctrl = Acc_math(idx_ctrl); 
Acc_str_ctrl = Acc_str(idx_ctrl);

figure,
subplot 211
bar([Acc_str_patn;Acc_math_patn]', 0.4), title('Acc-patn')
set(gca,'Xtick', 1:length(Acc_str_patn),'XtickLabel',1:length(Acc_str_patn));
xlabel('Sub'), ylabel('Acc (%)'),
set(gca,'FontSize',10,'XTickLabelRotation',90);
grid on
set(gca,'color','none');
legend({'Str','Math'})

subplot 212
bar([Acc_math_ctrl;Acc_str_ctrl]', 0.4), title('Acc-ctrl')
set(gca,'Xtick', 1:length(Acc_math_ctrl),'XtickLabel',1:length(Acc_math_ctrl));
xlabel('Sub'), ylabel('Acc (%)'),
set(gca,'FontSize',10,'XTickLabelRotation',90);
grid on
set(gca,'color','none');
legend({'Str','Math'})

%%
figure,
subplot 211
bar([mean(Acc_str_patn); mean(Acc_math_patn)]', 0.4), title('Acc-patn')
set(gca,'Xtick', 1:2,'XtickLabel',{'Str','Math'});
xlabel('Sub'), ylabel('Acc (%)'),
% set(gca,'FontSize',10,'XTickLabelRotation',90);
grid on
set(gca,'color','none');
legend({'Str','Math'})

subplot 212
bar([mean(Acc_math_ctrl); mean(Acc_str_ctrl)]', 0.4), title('Acc-ctrl')
set(gca,'Xtick', 1:2,'XtickLabel',{'Str','Math'});
xlabel('Sub'), ylabel('Acc (%)'),
% set(gca,'FontSize',10,'XTickLabelRotation',90);
grid on
set(gca,'color','none');
legend({'Str','Math'})

%%
cd('/data/MEG/Projects/ECP/behav')
figure,
bar([mean(Acc_str_patn); mean(Acc_math_patn)]', 0.4), title('Acc-patn')
hold on
errorbar(1:2,[mean(Acc_str_patn); mean(Acc_math_patn)]',[std(Acc_str_patn); std(Acc_math_patn)]')
set(gca,'Xtick', 1:2,'XtickLabel',{'Str','Math'});
xlabel('Sub'), ylabel('Acc (%)'),
% set(gca,'FontSize',10,'XTickLabelRotation',90);
grid on
set(gca,'color','none');
% legend({'Str';'Math'})

disp([mean(Acc_str_patn); mean(Acc_math_patn)])
disp([std(Acc_str_patn); std(Acc_math_patn)])

%% Response-time
clc
outdir = '/rcc/stor1/projects/ECP/MEG/';
outd_sub = rdir(fullfile(outdir,'MEG_work_ft','**', 'RT*.mat'));

rt = [];
for i=1:length(outd_sub)
    load(outd_sub(i).name);
    tkz = tokenize(outd_sub(i).name,'/');
    rt.all(i) = (m_RT_str_uniq+m_RT_math_uniq)/2; rt.str(i) = m_RT_str_uniq; rt.math(i) = m_RT_math_uniq;
    sub_name{i} = tkz{end-3};
%     sub_run{i} = tkz{end-1};
end

sub_name1 = unique(sub_name);

% find(rt.all) = nan;
% [row, col] = find(isnan(rt.all));
% rt.all(col)
% rt.math(col)
% rt.str(col)
% sub_name{col}

rt_runs = [];
for i=1:length(sub_name1)
    idx = find(strcmp(sub_name1{i},sub_name)==1);
    if length(idx) ==2
    rt_runs.all(i) = (rt.all(idx(1)) + rt.all(idx(2)))/2; 
    rt_runs.str(i) = (rt.str(idx(1)) + rt.str(idx(2)))/2; 
    rt_runs.math(i) = (rt.math(idx(1)) + rt.math(idx(2)))/2;
    tmp = sub_name(idx(1));
    sub_runs_name{i} = tmp{1};
    elseif length(idx) ==1
    rt_runs.all(i) = rt.all(idx); 
    rt_runs.str(i) = rt.str(idx); 
    rt_runs.math(i) = rt.math(idx);
    tmp = sub_name(idx);
    sub_runs_name{i} = tmp{1};
    disp(sub_name(idx))
    end
end

figure,
bar(rt_runs.all), L = length(rt_runs.all);
set(gca,'Xtick', 1:L,'XtickLabel',sub_runs_name);
set(gca,'FontSize',10,'XTickLabelRotation',90);
set(gcf, 'Position', [1000   100   1500   300]);
% idx_badperf = rt_runs.all < 50;
% disp(sub_runs_name(idx_badperf));
set(gca,'color','none');

%%
[~, ia,ib] = intersect(sub_ID,sub_runs_name);

sub_ID1 = sub_ID(ia);
RT_str = rt_runs.str(ib);
RT_math = rt_runs.math(ib);

RT_math_patn = RT_math(idx_patn); 
RT_str_patn = RT_str(idx_patn); 

RT_math_ctrl = RT_math(idx_ctrl); 
RT_str_ctrl = RT_str(idx_ctrl);

figure,
subplot 211
bar([RT_str_patn;RT_math_patn]', 0.4), title('RT-patn')
set(gca,'Xtick', 1:length(RT_str_patn),'XtickLabel',1:length(RT_str_patn));
xlabel('Sub'), ylabel('RT (sec)'),
set(gca,'FontSize',10,'XTickLabelRotation',90);
grid on
set(gca,'color','none');
legend({'Str','Math'})

subplot 212
bar([RT_math_ctrl;RT_str_ctrl]', 0.4), title('RT-ctrl')
set(gca,'Xtick', 1:length(RT_math_ctrl),'XtickLabel',1:length(RT_math_ctrl));
xlabel('Sub'), ylabel('RT (sec)'),
set(gca,'FontSize',10,'XTickLabelRotation',90);
grid on
set(gca,'color','none');
legend({'Str','Math'})

%%
figure,
subplot 211
bar([nanmean(RT_str_patn); nanmean(RT_math_patn)]', 0.4), title('RT-patn')
set(gca,'Xtick', 1:2,'XtickLabel',{'Str','Math'});
xlabel('Sub'), ylabel('RT (sec)'),
% set(gca,'FontSize',10,'XTickLabelRotation',90);
grid on
set(gca,'color','none');
legend({'Str','Math'})

subplot 212
bar([nanmean(RT_math_ctrl); nanmean(RT_str_ctrl)]', 0.4), title('RT-ctrl')
set(gca,'Xtick', 1:2,'XtickLabel',{'Str','Math'});
xlabel('Sub'), ylabel('RT (sec)'),
% set(gca,'FontSize',10,'XTickLabelRotation',90);
grid on
set(gca,'color','none');
legend({'Str','Math'})

%%





