%% ECP Semantic decision task dataset, Medical College of Wisconsin

% Script: Process (extract virtual sensrors at voxel and roi AAL atlas levels)
% Project: ECP_SD
% Writtern by: Vahab Youssof Zadeh
% Update: 11/16/2022

clear; clc, close('all'); warning off

%% Paths
addpath('./run')
Run_setpath
addpath('./data')
addpath('./run')
addpath(genpath('./functions'))

% clear; clc, close('all'); warning off
% 
% %%
% cd '/data/MEG/Vahab/Github/MCW_MEGlab/FT';
% restoredefaultpath
% cd_org = cd;
% addpath(genpath(cd_org));
% % rmpath('./Failedattemps');
% 
% indir = '/data/MEG/Research/Aqil_Izadi/process/MNE_analysis/dSPM';
indir = '/data/MEG/Research/ECP/Semantic_Decision/BS_database/';
% outdir = '/data/MEG/Research/Aqil_Izadi/process/MNE_analysis/analysis';
outdir = '/data/MEG/Research/ECP/Semantic_Decision/process/process_lcmv';
% 
% %- Adding path
% cfg_init = [];
% cfg_init.path_tools = '/data/MEG/Vahab/Github/MCW_MEGlab/tools';
% [allpath, atlas] = vy_init(cfg_init);
% 
% addpath('/data/MEG/MCW pipelines/functions');

%% Reading LCMV data, event 2, Symbols
clear subj_all
pardir = 'Parcellation/HCPMMP1/Mean/parc/HCPMMP1/event_2/sub_by_sub/';
d = rdir(fullfile (indir, pardir,'/**/*.csv'));
clear val subj run
for i=1:length(d)
    disp(num2str(i))
    [pathstr, name] = fileparts(d(i).name);
    datafolder{i} = pathstr;
    datafile{i} = d(i).name;
    tkl = tokenize(pathstr,'/');
    name = tkl{end};
    Index = strfind(name, 'ec'); Index1 = strfind(name, 'run');
    subj_all{i} = [num2str(i), ': ', name(Index(1):Index(1)+5), '_run',name(Index1(1)+3)];
    subj{i} = name(Index(1):Index(1)+5);
    run{i} = name(Index1(1)+3);
    val_2{i} = readtable(datafile{i});
end

pardir = 'Parcellation/HCPMMP1/Mean/parc/HCPMMP1/event_3/sub_by_sub/';
d = rdir(fullfile (indir, pardir,'/**/*.csv'));
clear val
for i=1:length(d)
    disp(num2str(i))
    [pathstr, name] = fileparts(d(i).name);
    datafolder{i} = pathstr;
    datafile{i} = d(i).name;
    tkl = tokenize(pathstr,'/');
    name = tkl{end};
    Index = strfind(name, 'ec'); Index1 = strfind(name, 'run');
    subj_all{i} = [num2str(i), ': ', name(Index(1):Index(1)+5), '_run',name(Index1(1)+3)];
    subj{i} = name(Index(1):Index(1)+5);
    run{i} = name(Index1(1)+3);
    val_3{i} = readtable(datafile{i});
end

%%
cd(outdir)
load('dSPM3.mat'), val_3 = val;
load('dSPM2.mat'), val_2 = val;

%%
src_fname = '/data/MEG/Vahab/Github/MCW_MEGlab/FT/Data_file/template/mne_fsaverage/fsaverage/bem/fsaverage-ico-5-src.fif';
src = ft_read_headshape(src_fname, 'format', 'mne_source');

src_S = [];
src_s.pnt = src.pos;
src_s.tri = src.tri;

figure,
ft_plot_mesh(src_s, 'faecolor', 'brain', 'facealpha', 1);

%%
num_sub = length(val_2);
for i = 1:num_sub
    val_symb2(i,:,:) = val_2{i};
    val_anim3(i,:,:) = val_3{i};
end

%%
disp('1: dSPM animal')
disp('2: dSPM symbol')
disp('3: dSPM animal-symbol')
d_Sel = input(':');

switch d_Sel
    case 1
        d_in = val_anim3; lbl = 'Animal, 3';
    case 2
        d_in = val_symb2;  lbl  = 'Symbol, 2';
    case 3
        d_in = val_anim3 - val_symb2; lbl  = 'Animal - Symbol, 3-vs-2';
end

%%
d_in_sub = (d_in(1:2:end,:,:) + d_in(2:2:end,:,:))./2;
size(d_in_sub)

% subj_all(1:2:end)
% subj_all(2:2:end)

%% plot mean
grand_val  = squeeze(mean(d_in,1));

grand_val_3  = squeeze(mean(val_anim3,1));
grand_val_2  = squeeze(mean(val_symb2,1));

%
figure,
ft_plot_mesh(src_s, 'vertexcolor', mean(grand_val,2));
colorbar, view(-90, 0); title(lbl)

figure
cfg = [];
cfg.view = [-180,-90;0,90;-90,0; 90,0];
cfg.position = [800   800   1000   300];
cfg.color = (viridis(256));
mcw_surfplot(cfg,src_s,mean(grand_val,2));

Time = linspace(0,1,size(grand_val,2));

%% Plot mean values (over time)
% in time
close all
thre = 0;
figure,
cfg = [];
% cfg.view = [-180,-90;0,90;-90,0; 90,0];
cfg.view = [-90,0; 90,0];
cfg.color = (viridis(256));
% cfg.h = ['L','R'];
cfg.position = [800   800   1000   300];
for i=1:20:size(grand_val,2)
    tmp  = grand_val(:,i);
    tmp(tmp < thre.*max(tmp(:))) = 0;
%     tmp = tmp.^0.25;
    mcw_surfplot(cfg,src_s, tmp);
    title([num2str(Time(i)), ' s'])
    pause(1)
end

%% Plot mean values (over time intervals)
w1 = 0; l = 0.1; ov = l.*1; j=1; wi=[];
while w1+l <= 1
    wi(j,:) = [w1, w1+l]; j=j+1; w1 = w1 + ov;
end
disp(wi)
pause,

close all
thre = 0;
figure,
cfg = [];
% cfg.view = [-180,-90;0,90;-90,0; 90,0];
cfg.view = [-90,0; 90,0];
cfg.color = (viridis(256));
% cfg.h = ['L','R'];
cfg.position = [800   800   1000   300];
for i=1:length(wi)
    
    [~, idx_toi1] = min(abs(Time - wi(i,1))); 
    [val, idx_toi2] = min(abs(Time - wi(i,2)));
    
    tmp  = mean(grand_val(:,idx_toi1:idx_toi2),2);
    tmp(tmp < thre.*max(tmp(:))) = 0;
%     tmp = tmp.^0.25;
    mcw_surfplot(cfg,src_s, tmp);
    title([num2str(wi(i,1)), '-', num2str(wi(i,2)), ' sec'])
    pause(1)
end

%% Plot mean values (single time intervals)
thre = 0;
toi_pst = input('enter time of interests, e.g.,[0,1]: ');
[~, idx_toi1] = min(abs(Time - toi_pst(1))); [val, idx_toi2] = min(abs(Time - toi_pst(2)));
tmp  = mean(grand_val(:,idx_toi1:idx_toi2),2);
tmp(abs(tmp) < thre.*max(tmp(:))) = 0;
%     tmp = tmp.^0.25;
cfg = [];
cfg.view = [-180,-90;0,90;-90,0; 90,0];
cfg.color = (viridis(256));
cfg.position = [800   800   1000   300];
mcw_surfplot(cfg,src_s, tmp);

% spmpath = '/data/MEG/Vahab/Github/MCW_MEGlab/tools/SPM/spm12_2021/spm12';
% addpath(spmpath)
% spm_get_defaults
% conn_path = '/data/MEG/Vahab/Github/MCW_MEGlab/tools/Conn/conn';
% addpath(conn_path);

% conn_mesh_display(src_s, '', '', '', tmp);
% title(par)

%% Plot stats (over time intervals)
w1 = 0; l = 0.1; ov = l.*1; j=1; wi=[];
while w1+l <= 1
    wi(j,:) = [w1, w1+l]; j=j+1; w1 = w1 + ov;
end
disp(wi)
pause,

val_in = grand_val;

close all
thre = 0;
figure,
cfg = [];
cfg.view = [-180,-90;0,90;-90,0; 90,0];
% cfg.view = [-90,0; 90,0];
% cfg.color = flipud(viridis(256));
cfg.color = (viridis(256));
% cfg.h = ['L','R'];
cfg.position = [800   800   1000   300];
for i=1:length(wi)
    
    [~, idx_toi1] = min(abs(Time - wi(i,1))); 
    [val, idx_toi2] = min(abs(Time - wi(i,2)));
    
    tmp  = mean(val_in(:,idx_toi1:idx_toi2),2);
%     tmp(tmp < thre.*max(tmp(:))) = 0;
%     tmp = tmp.^0.25;
    mcw_surfplot(cfg,src_s, tmp);
    title([num2str(wi(i,1)), '-', num2str(wi(i,2)), ' sec'])
    pause(1)
end

%% Correlation (against mean)
num_sub = length(val_2);
clear m_val_2 m_val_3 val_sel2 val_sel3
for j = 1:201 
    clear stat_out
    disp([num2str(j), '/', num2str(size(val_2{1},2))])
    for i = 1:num_sub
        val_sel2(i,:) = val_2{i}(:,j);
        val_sel3(i,:) = val_3{i}(:,j);
    end
    m_val_2(j,:) = mean(val_sel2,1);
    m_val_3(j,:) = mean(val_sel3,1);
    
%     size(m_val_2)
%     size(val_sel2)
    
    for i=1:size(val_sel2,1)
        crr_2(j,i) = corr(m_val_2(j,:)',val_sel2(i,:)', 'type','Spearman');
        crr_3(j,i) = corr(m_val_3(j,:)',val_sel3(i,:)', 'type','Spearman');
    end   
end

%%
close all
%- RUN-LEVEL Corr
% figure, 
% subplot 211, imagesc(crr_2'), title('SYMB'), ylabel('RUN')
% subplot 212, imagesc(crr_3'), title('ANIM'), ylabel('RUN')

%- SUBJECT-LEVEL Corr
sub_crr_2 = (crr_2(:,1:2:end) + crr_2(:,2:2:end))./2; 
figure, 
subplot 211
imagesc(sub_crr_2'), title('SYMB'), ylabel('SUB')
subplot 212
sub_crr_3 = (crr_3(:,1:2:end) + crr_3(:,2:2:end))./2; 
imagesc(sub_crr_3'), title('ANIM'), ylabel('SUB')

%%
% figure, bar(nanmean(sub_crr_2,1))
% figure, bar(nanmean(sub_crr_3,1))
corr_23_sub = (nanmean(sub_crr_2,1)  + nanmean(sub_crr_3,1))./2;

figure,
subplot 211
bar(corr_23_sub); xlabel('subj'); ylabel('mean corr'); title('ANIM & SYMB'); 
set(gca,'color','none');
set(gca,'Xtick', 1:length(corr_23_sub))
set(gca,'FontSize',10,'XTickLabelRotation',90);

corr_23_time = (nanmean(sub_crr_2,2)  + nanmean(sub_crr_3,2))./2;
% figure, 
subplot 212
plot(Time,corr_23_time); xlabel('time (sec)'); ylabel('mean corr'); title('ANIM & SYMB');  
set(gca,'color','none');


%%
[I, SUB_best] = sort(corr_23_sub, 'descend');

% size(corr_23_sub)
% size(d_in_sub)
sub_val = mean(d_in_sub(SUB_best(1),:,:),3)';
% sub_val = mean(d_in_sub(14,:,:),3)';
% size(sub_val)

% figure,
% ft_plot_mesh(src_s, 'vertexcolor',sub_val');

figure
cfg = [];
cfg.view = [-180,-90;0,90;-90,0; 90,0];
cfg.position = [800   800   1000   300];
cfg.color = (viridis(256));
mcw_surfplot(cfg,src_s,sub_val);

%% Response-time
clc,
rtdir = '/data/MEG/Research/Aqil_Izadi/process';
d = rdir([rtdir,'/**/*RT.mat']);

RT_all = [];
for j=1:length(d)    
    RT_sel = load(d(j).name);
    RT_all.sub{j} = RT_sel.subj;
    RT_all.run{j} = RT_sel.run;
    RT_all.rt_time.animal(j) = mean(RT_sel.rt_time.animal);
    RT_all.rt_time.symbol(j) = mean(RT_sel.rt_time.symbol);
    RT_all.rt_time.both(j) = mean(RT_sel.rt_time.both);   
end

%%
cd(fullfile(indir, '/event_3/dspm_stc_morph_dict_event 3_csv_data'))
datafile1 = [];
d = rdir('./*_stc*.csv');
for i=1:length(d)
    disp(num2str(i))
    [pathstr, name] = fileparts(d(i).name);
    datafolder{i} = pathstr;
    datafile{i} = d(i).name;
    Index = strfind(name, 'ec'); Index1 = strfind(name, 'run');
    subj_ID{i} = ['EC', name(Index(1)+2:Index(1)+5)];
end

%%
[C,IA,IB] =  intersect(RT_all.sub,subj_ID');

%%
idx_val = [];
for i=1:length(C)
    idx = find(strcmp(C{i},subj_ID)==1);
    idx_val = [idx_val,idx];
end

idx_rt = [];
for i=1:length(C)
    idx = find(strcmp(C{i},RT_all.sub)==1);
    idx_rt = [idx_rt,idx];
end

%%
rt_3 = RT_all.rt_time.animal(idx_rt);
rt_2 = RT_all.rt_time.symbol(idx_rt);

crr_2_sel = crr_2(:,idx_val);
crr_3_sel = crr_3(:,idx_val);

close all
sub_crr_sel_2 = (crr_2_sel(:,1:2:end) + crr_2_sel(:,2:2:end))./2; figure, imagesc(sub_crr_sel_2'), title('2')
sub_crr_sel_3 = (crr_3_sel(:,1:2:end) + crr_3_sel(:,2:2:end))./2; figure, imagesc(sub_crr_sel_3'), title('3')

sub_rt_2 = (rt_2(:,1:2:end) + rt_2(:,2:2:end))./2; figure, imagesc(sub_rt_2'), title('2')
sub_rt_3 = (rt_3(:,1:2:end) + rt_3(:,2:2:end))./2; figure, imagesc(sub_rt_3'), title('3')

x = nanmean(sub_crr_sel_2,1); y = sub_rt_2;
figure, plot(x,y,'*'), hold on; P = polyfit(x,y,1); disp(P(1)),[r,p] = corrcoef(x,y)
x_min = min(x(:)); x_max = max(x(:)); d_min = polyval(P,0); d_max = polyval(P,1); plot([0 1],[d_min d_max],'k--')

x = nanmean(sub_crr_sel_3,1); y = sub_rt_3;
figure, plot(x,y,'*'), hold on; P = polyfit(x,y,1); disp(P(1)),[r,p] = corrcoef(x,y)
x_min = min(x(:)); x_max = max(x(:)); d_min = polyval(P,0); d_max = polyval(P,1); plot([0 1],[d_min d_max],'k--')

x = nanmean(crr_2_sel,1); y = rt_2;
figure, plot(x,y,'*'), hold on; P = polyfit(x,y,1); disp(P(1)),[r,p] = corrcoef(x,y)
x_min = min(x(:)); x_max = max(x(:)); d_min = polyval(P,0); d_max = polyval(P,1); plot([0 1],[d_min d_max],'k--')

x = nanmean(crr_3_sel,1); y = rt_3;
figure, plot(x,y,'*'), hold on; P = polyfit(x,y,1); disp(P(1)),[r,p] = corrcoef(x,y)
x_min = min(x(:)); x_max = max(x(:)); d_min = polyval(P,0); d_max = polyval(P,1); plot([0 1],[d_min d_max],'k--')

x = (nanmean(crr_2_sel,1) + nanmean(crr_3_sel,1))./2; y = (rt_3 + rt_2)./2;
figure, plot(x,y,'*'), hold on; P = polyfit(x,y,1); disp(P(1)),[r,p] = corrcoef(x,y)
x_min = min(x(:)); x_max = max(x(:)); d_min = polyval(P,0); d_max = polyval(P,1); plot([0 1],[d_min d_max],'k--')

%%
% % pyversion
% % 
% % x = py.numpy.linespace(0,10,101)
% % 
% % 
% % filename = '/data/MEG/Research/Aqil_Izadi/process/MNE_analysis/Processed MPVA/Second Analysis/avr_stc_morph.pkl';
% % fid=py.open(filename,'rb');
% % data=py.pickle.load(fid);
% % 
% % 
% % %%
% % addpath('/data/MEG/Research/Aqil_Izadi/Scripts/functions')
% % 
% % [a] = loadpickle(filename)
% % 
% % %%
% % filename = '/data/MEG/Research/Aqil_Izadi/process/MNE_analysis/Processed MPVA/Second Analysis/avr_stc_morph.pkl';
% % fid = py.open(filename,'rb');
% % data = py.pickle.load(fid);
% 
% clear,
% clc
% %%
% cd('/data/MEG/Research/Aqil_Izadi/process/MNE_analysis/Processed MPVA/Second Analysis')
% 
% %% read .cvs files (tgm)
% filename = 'grand tgm average.csv';
% M = csvread(filename);
% 
% figure,
% imagesc(M)
% 
% 
% %% read .cvs files (tgm)
% filename = 'raws_with_event_exception.csv';
% M1 = csvread(filename);
% 
% figure,
% imagesc(M)
% 
% 
% %% First analysis
% cd('/data/MEG/Research/Aqil_Izadi/process/MNE_analysis/Processed MPVA/First Analysis/MVPA INVs')
% cd('/data/MEG/Research/Aqil_Izadi/process/MNE_analysis/Processed MPVA/First Analysis/MVPA STCs')
% cd('/data/MEG/Research/Aqil_Izadi/process/MNE_analysis/Processed MPVA/First Analysis/MVPA EVs')
% 
% hdr = ft_read_header('ec1105_SD_run1_raw-inv.fif'); % your fif-filename
% % hdr = ft_read_header('ec2086_SD_run2_raw-rh.stc'); % your fif-filename
% hdr = ft_read_header('ec1091_SD_run2_raw-ave.fif'); % your fif-filename
% 
% 
% cfg = [];
% cfg.dataset = 'ec1091_SD_run2_raw-ave.fif'
% data = ft_preprocessing(cfg);


