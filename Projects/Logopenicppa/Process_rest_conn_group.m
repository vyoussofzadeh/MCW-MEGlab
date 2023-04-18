%% Logopenicppa dataset, Medical College of Wisconsin

% Script: BS Process (seed-based conn analysis)
% Project: Logopenicppa_rest
% Writtern by: Vahab YoussofZadeh
% Update: 04/12/2023

clear; clc, close('all'); warning off

%% Analysis flags
flag = [];
flag.plottriggers = 1;
flag.anatomy_check = 1;
flag.dicsanalysis = 0;
flag.connanalysis = 1;

%% Paths
% addpath('./run')
% Run_setpath
addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/Logopenicppa/functions')
[atlas, all_path] = do_setpath();

%%
indir = '/data/MEG/Research/logopenicppa';
outdir = fullfile(indir,'ft_process');
rawdatadir = fullfile(indir,'raw');

%%
clc
ft_progress('init', 'text',     'please wait ...');
d  = rdir(fullfile(outdir,'/**/Rest/**/seed*.mat'));
clear Sub Run Session
for i=1:length(d)
    ft_progress(i/length(d), 'Processing network values %d from %d', i, length(d));
    tkz = tokenize(d(i).name,'/');
    conn = load(d(i).name);
    network_atlas{i} = conn.net_conn_seed.network_atlas;
    network{i} = conn.net_conn_seed.network;
    source_conn{i} = conn.net_conn_seed.source_conn;
    pause(0.1);
    Sub{i} = tkz{end-3};
    Session{i} = tkz{end-1};
    run{i} = tkz{end}(end-4);
end
ft_progress('close');

%%
for i=1:length(d)
    cfg = [];
    cfg.subj = [Sub{i}, '-', Session{i}, '-', run{i}];
    cfg.mask = 'tissue';
    cfg.thre = 0;
    cfg.savepath = '/data/MEG/Research/logopenicppa/results/SeedConn'; %'groupave';
    cfg.colorbar = 2;
    cfg.saveflag = 1;
    % cfg.colormap = colormap(flipud(pink));
    cfg.colormap = [];
    cfg.surfinflated   = 'surface_inflated_both.mat';
%     cfg.views = [90 0;0,90;-90 0;-180,-90];
    cfg.views = [-90,0; 90,0];
    cfg.tit = [Sub{i}, '-', Session{i}, '-', run{i}];
    do_mapvis(cfg, network_atlas{i});
    pause(0.1),
    close all
end

%%
% mnetwork_atlas = network_atlas{1};
close all
k=1; k1=1;
clear conn_sub_all
for i=[2,4,5,6]
    idx = find(strcmp(Sub, ['32037_00', num2str(i)])==1);
    disp(idx)
    %     mnetwork_atlas = network_atlas{1};
    k=1;
    clear conn_sub
    for j=idx
        %         net_sub = network_atlas {j};
        disp(j)
        conn_sub(:,:,:,k) = source_conn{j}.wpli_debiasedspctrm;
        k=k+1;
    end
    conn_sub_all{k1} = conn_sub;
    k1=k1+1;
end

mm = mean((conn_sub_all{1}),4);

foi = [1,30];
[~, idx_mn] = min(abs(foi(1) - source_conn{1}.freq));
[~, idx_mx] = min(abs(foi(2) - source_conn{1}.freq));

% aedge =  squeeze(source_conn1.(par)(:,:,1));
aedge =  mean(mm(:,:,idx_mn:idx_mx),3);
aedge(nonIdenticalIndices,:) = 0;

%%
% mnetwork_atlas = network_atlas{1};
close all
sub_id = [2,4,5,6];
k=1; k1=1;
network_atlas_sub = [];
clear network_atlas_sub
for i = sub_id
    idx = find(strcmp(Sub, ['32037_00', num2str(i)])==1);
    disp(idx)
    %     mnetwork_atlas = network_atlas{1};
    k=1;
    clear conn_sub
    for j=idx
        %         net_sub = network_atlas {j};
        disp(j)
        network_atlas_sub(:,:,:,k) = network_atlas{j}.tissue;
        k=k+1;
    end
    network_atlas_sub_all{k1} = network_atlas_sub;
    k1=k1+1;
end

sel = 1;
mm = mean((network_atlas_sub_all{sel}),4);
mm(mm < 0.7.*max(mm(:))) = 0;
mm = round(mm);

% mm = (squeeze(network_atlas_sub_all{1}(:,:,:,4))); mm = round(mm);

mnetwork_atlas = network_atlas{1};
mnetwork_atlas.tissue = mm;

cfg = [];
cfg.subj = ['32037_00', num2str(sub_id(sel))];
cfg.mask = 'tissue';
cfg.thre = 0;
cfg.savepath = 'groupave';
cfg.colorbar = 0;
cfg.saveflag = 0;
% cfg.colormap = colormap(flipud(pink));
cfg.colormap = [];
cfg.surfinflated   = 'surface_inflated_both.mat';
%     cfg.views = [90 0;0,90;-90 0;-180,-90];
cfg.views = [-90,0; 90,0];
cfg.tit = ['32037_00', num2str(sub_id(sel))];
do_mapvis(cfg, mnetwork_atlas);

%%
% mnet = [];
% k=1;
% for j=idx
% %         net_sub = network_atlas {j};
%         disp(j)
%         tmp = network{j}.eigenvector_cent;
%         mnet(:,:,k) =  tmp;
%         k=k+1;
%     end
% %     tmp = mnet./ length(idx);
% %     
% %     tmp = (tmp - min(tmp(:))) ./ (max(tmp(:)) - min(tmp(:)));
%     
% %     mnetwork_atlas.tissue = mnet;
%     
%     cfg = [];
%     cfg.subj = ['32037_00', num2str(i)];
%     cfg.mask = 'tissue';
%     cfg.thre = 0;
%     cfg.savepath = 'groupave';
%     cfg.colorbar = 2;
%     cfg.saveflag = 0;
%     % cfg.colormap = colormap(flipud(pink));
%     cfg.colormap = [];
%     cfg.surfinflated   = 'surface_inflated_both.mat';
%     %     cfg.views = [90 0;0,90;-90 0;-180,-90];
%     cfg.views = [-90,0; 90,0];
%     cfg.tit = ['32037_00', num2str(i)];
%     do_mapvis(cfg, mnetwork_atlas);
% end
% 
% %%
% close all
% savedir = '/group/bgross/work/CIDMEG/analysis/process/pilot1/process_100Hz';
% ddd = rdir(fullfile(savedir, '/**/wPLI*.mat'));
% 
% clear ID conn_val
% for i=1:length(ddd)
%     tkz = tokenize(ddd(i).name,'/');
%     conn = load(ddd(i).name);
%     conn_val(:,:,:,i) = conn.source_conn1.wpli_debiasedspctrm;
%     ID{i} = tkz{end}(1:end-4);
% end
% mconn_val = nanmean(conn_val,4);
% 
% source_conn1 = conn.source_conn1;
% source_conn1.wpli_debiasedspctrm = mconn_val;
% par = conn.par;
% vs_roi1 = conn.vs_roi1;
% 
% cd('/group/bgross/work/CIDMEG/analysis/process/pilot1/process_100Hz/group')
% Run_networkanalysis %- network analysis (network degrees)
% Run_vis_netconn_analysis %- mapping
% 
% 
% %% Anatomy
% subj = 'pilot1';
% run = 1;
% datafilename = ['Task_run', num2str(run), '_raw.fif'];
% datafile = fullfile(indir, datafilename);
% subdir = fullfile(outdir, subj,'preprocess_100Hz'); % output dir
% load(fullfile(subdir,['data_run', num2str(run), '_cond',num2str(1), '.mat']))
% 
% mridir = '/group/bgross/work/CIDMEG/analysis/anatomy/Sub001/T1';
% mrifile = 'sub-CIDFMRI001_ses-1_acq-mprage_T1w.nii.gz';
% Run_anatomy