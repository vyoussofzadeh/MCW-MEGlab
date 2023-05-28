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
conn_method = 'amplcorr'; par = [conn_method, 'spctrm'];
conn_method = 'wpli_debiased'; par = [conn_method, 'spctrm'];

%%
clc
ft_progress('init', 'text',     'please wait ...');
% d  = rdir(fullfile(outdir,'/**/Rest/**/seed*1-30Hz.mat'));
d  = rdir(fullfile(outdir,['/**/Rest/**/seed*',conn_method, '.mat']));
clear Sub Run Session source_conn network network_atlas Seed_idx
for i=1:length(d)
    ft_progress(i/length(d), 'Processing network values %d from %d', i, length(d));
    disp(d(i).name)
    tkz = tokenize(d(i).name,'/');
    conn_comput = load(d(i).name);
    network_atlas{i} = conn_comput.net_conn_seed.network_atlas;
    network{i} = conn_comput.net_conn_seed.network;
    source_conn{i} = conn_comput.net_conn_seed.source_conn;
    pause(0.1);
    Sub{i} = tkz{end-3}; Session{i} = tkz{end-1};
    tkz1 = tokenize(tkz{end},'_');
    run{i} = tkz1{end-1};
    Seed_idx{i} = conn_comput.net_conn_seed.seed_idx;
end
ft_progress('close');

%%
datainfo = [];
datainfo.Seed_idx = Seed_idx;
datainfo.Sub = Sub;
datainfo.Session = Session;
datainfo.run = run;
datainfo.UQSeed_idx = unique(cell2mat(datainfo.Seed_idx));
datainfo.UQSub = unique((datainfo.Sub));

%%
for i=1:length(d)
    cfg = [];
    cfg.subj = [Sub{i}, '-', Session{i}, '-', run{i}];
    cfg.mask = 'tissue';
    cfg.thre = 0;
    cfg.savepath = '/data/MEG/Research/logopenicppa/results/1_30/SeedConn'; %'groupave';
    cfg.colorbar = 2;
    cfg.saveflag = 0;
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
        conn_sub(:,:,:,k) = source_conn{j}.(par);
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
% aedge(nonIdenticalIndices,:) = 0;

%% Map seed
close all
unq_seed_idx = unique(cell2mat(Seed_idx));
for j=1:1%length(unq_seed_idx)
    seed_map = atlas;
    seed_map.tissue = zeros(size(atlas.tissue));
    for i = 1:116
        idx = atlas.tissue == unq_seed_idx(j);
        seed_map.tissue(idx) = 1;
    end
    
    cfg = [];
    cfg.subj = ['seed:', atlas.tissuelabel{unq_seed_idx(j)}];
    cfg.mask = 'tissue';
    cfg.thre = 0;
    cfg.savepath = '/data/MEG/Research/logopenicppa/results/1_30/Seed_roi';
    cfg.colorbar = 0;
    cfg.saveflag = 1;
    cfg.colormap = [];
    cfg.surfinflated   = 'surface_inflated_both.mat';
    % cfg.views = [90 0;0,90;-90 0;-180,-90];
    cfg.views = [-90,0; 90,0];
    cfg.tit = ['seed:', atlas.tissuelabel{unq_seed_idx(j)}];
    do_mapvis(cfg, seed_map);    
end

%% Subj - group
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

sel = 2;
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

%% Subj, session - group
close all
sub_id = [2,4,5,6]; k1=1;

clear network_atlas_sub network_atlas_sub_all
for i = sub_id
    idx = find(strcmp(Sub, ['32037_00', num2str(i)])==1);
    disp(idx)
    %     mnetwork_atlas = network_atlas{1};
    
    suq = unique(Session(idx));
    clear conn_sub network_atlas_sub
    for kk=1:length(suq)
        idx2 = find(strcmp(suq(kk), Session(idx))==1);
        k=1;
        for j=idx(idx2)
            disp(j)
            network_atlas_sub(:,:,:,k,kk) = network_atlas{j}.tissue;
            k=k+1;
        end
    end
    network_atlas_sub_all{k1} = network_atlas_sub;
    k1=k1+1;
end

size(network_atlas_sub_all{2})

for sel = 1:1%length(sub_id)
    for kk=1:size(network_atlas_sub_all{sel},5)
        
        mm = network_atlas_sub_all{sel}(:,:,:,:,kk);
        size(mm)
        mm = mean(mm,4);
        mm(mm < 0.7.*max(mm(:))) = 0;
        mm = round(mm);
        
        mnetwork_atlas = network_atlas{1};
        mnetwork_atlas.tissue = mm;
        
        cfg = [];
        cfg.subj = ['32037_00', num2str(sub_id(sel)), ' session: ', num2str(kk)];
        cfg.mask = 'tissue';
        cfg.thre = 0;
        cfg.savepath = '/data/MEG/Research/logopenicppa/results/1_30/Seed_conn_session';
        cfg.colorbar = 0;
        cfg.saveflag = 1;
        % cfg.colormap = colormap(flipud(pink));
        cfg.colormap = [];
        cfg.surfinflated   = 'surface_inflated_both.mat';
        %     cfg.views = [90 0;0,90;-90 0;-180,-90];
        cfg.views = [-90,0; 90,0];
        cfg.tit = ['32037_00', num2str(sub_id(sel)), ' session: ', num2str(kk)];
        do_mapvis(cfg, mnetwork_atlas);
        pause(0.1),
        close all
    end
end

%%
addpath('/data/MEG/Vahab/Github/MCW_MEGlab/tools/Conn/conn');
h = get(0, 'Children');
if isempty(findobj(h,'tag','CONN functional connectivity toolbox'))
    conn;
end
path_tools = '/data/MEG/Vahab/Github/MCW_MEGlab/tools';
spm_path = fullfile(path_tools,'SPM/spm12');
addpath(genpath(spm_path))
spm_get_defaults

%% Plot conn - all runs, sessions
d  = rdir(fullfile(outdir,'/**/Rest/**/seed*1-30Hz.mat'));
for i=1:length(d)
    conn_comput = load(d(i).name);
    foi = conn_comput.net_conn_seed.foi;
    source_conn = conn_comput.net_conn_seed.source_conn;
%     par = 'wpli_debiasedspctrm';
    seed_idx = conn_comput.net_conn_seed.seed_idx;
    
    [~, idx_mn] = min(abs(foi(1) - source_conn.freq));
    [~, idx_mx] = min(abs(foi(2) - source_conn.freq));
    aedge =  mean(source_conn.(par)(:,:,idx_mn:idx_mx),3);
    
    tkz = tokenize(d(i).name,'/');
    Sub{i} = tkz{end-3}; Session{i} = tkz{end-1}; run{i} = tkz{end}(end-11);
    
    aedge1 = zeros(size(aedge));
    aedge1(seed_idx,:) = aedge(seed_idx,:);
    aedge1(:,seed_idx) = aedge(:,seed_idx);
    close all
    conn_mesh_display_vy('', '', '', Coordinate, aedge1, 0.2);
    disp(d(i).name)
    disp([Sub{i}, '_', Session{i}, '_', run{i},'.jpg'])
    savepath = '/data/MEG/Research/logopenicppa/results/1_30/Conn_maps';
    savefig = fullfile(savepath, [Sub{i}, '_', Session{i}, '_', run{i},'.jpg']);
    conn_print(savefig,...
        '-nogui',...
        '-row',...
        get(findobj(gcf,'label','Left view (both hem)'),'callback'),...
        get(findobj(gcf,'label','Anterior view'),'callback'),...
        get(findobj(gcf,'label','Superior view'),'callback'));
    pause(0.1),
    close all
end

%% Plot conn - mean-sessions
%- Subj, session - group
close all
sub_id = [2,4,5,6]; k1=1;
% sub_id = [6]; k1=1;

clear conn_comput_session_sub conn_comput_session_sub_all Seed_idx_session
for i = sub_id
    idx = find(strcmp(Sub, ['32037_00', num2str(i)])==1);
    disp(idx)
    %     mnetwork_atlas = network_atlas{1};
    
    suq = unique(Session(idx));
    clear conn_sub network_atlas_sub
    for kk=1:length(suq)
        idx2 = find(strcmp(suq(kk), Session(idx))==1);
        k=1;
        for j=idx(idx2)
            disp(j)
            conn_comput_session_sub(:,:,:,k,kk) = source_conn{j}.(par);
            Seed_idx_session(k,kk) = datainfo.Seed_idx{j};
            k=k+1;
        end
    end
    conn_comput_session_sub_all{k1} = conn_comput_session_sub;
    Seed_idx_session_all{k1} = unique(Seed_idx_session);
    k1=k1+1;
end
% size(conn_comput_session_sub_all{2})
% Seed_idx_session_all{1}

foi = [1,29];
for subsel = 1:length(datainfo.UQSub)
    for kk=1:size(conn_comput_session_sub_all{subsel},5)        
        aedge = conn_comput_session_sub_all{subsel}(:,:,:,:,kk);
        size(aedge)
        aedge = mean(aedge,4);
        [~, idx_mn] = min(abs(foi(1) - source_conn{1}.freq));
        [~, idx_mx] = min(abs(foi(2) - source_conn{1}.freq));
        aedge = abs(aedge);
        aedge =  mean(aedge(:,:,idx_mn:idx_mx),3);      
        seed_idx = Seed_idx_session_all{subsel};
        saveid = ['32037_00', num2str(sub_id(subsel)), ' session: ', num2str(kk)];        
        aedge1 = zeros(size(aedge));
        aedge1(seed_idx,:) = aedge(seed_idx,:);
        aedge1(:,seed_idx) = aedge(:,seed_idx);
        close all
%         conn_mesh_display_vy('', '', '', Coordinate, aedge1, 0.2);
        savepath = '/data/MEG/Research/logopenicppa/results/1_30/Seed_conn_maps';
        savefig = fullfile(savepath, [saveid,'.jpg']);
%         save(fullfile(savepath, [saveid,'.mat']), 'aedge1');
        save(fullfile(savepath, [saveid,'_', conn_method, '.mat']), 'aedge1');
%         conn_print(savefig,...
%             '-nogui',...
%             '-row',...
%             get(findobj(gcf,'label','Left view (both hem)'),'callback'),...
%             get(findobj(gcf,'label','Anterior view'),'callback'),...
%             get(findobj(gcf,'label','Superior view'),'callback'));
%         pause (0.1),
%         close all
    end
end

% Progression analysis
% d = rdir(fullfile(savepath, '*.mat'));
d = rdir(fullfile(savepath, ['*', conn_method, '*.mat']));


mm = [];
conn_mat_all = [];
for i=1:length(d)
    conn_mat = load(d(i).name);
    disp(d(i).name)
    conn_mat_all(i,:,:) = conn_mat.aedge1;
    mm(i) = nanmean(nanmean(conn_mat.aedge1(:)));
end

figure, bar(mm)

%%
temp = squeeze(mean(conn_mat_all(1:4,:,:),1)); size(temp)
conn_diff = squeeze(conn_mat_all(2,:,:) - conn_mat_all(1,:,:)); size(conn_diff)
conn_diff = squeeze(conn_mat_all(4,:,:) - conn_mat_all(3,:,:)); size(conn_diff)
conn_diff = squeeze(conn_mat_all(4,:,:) - temp); size(conn_diff)

conn_mesh_display_vy('', '', '', Coordinate, conn_diff, 0.2);

%%
clc, nanmean(nanmean(conn_diff(:)))

%% Plot conn - mean-sessions - whole brain
% par = 'wpli_debiasedspctrm';
%- Subj, session - group
close all
sub_id = [2,4,5,6]; k1=1;

clear conn_comput_session_sub conn_comput_session_sub_all Seed_idx_session
for i = sub_id
    idx = find(strcmp(Sub, ['32037_00', num2str(i)])==1);
    disp(idx)
    %     mnetwork_atlas = network_atlas{1};
    
    suq = unique(Session(idx));
    clear conn_sub network_atlas_sub
    for kk=1:length(suq)
        idx2 = find(strcmp(suq(kk), Session(idx))==1);
        k=1;
        for j=idx(idx2)
            disp(j)
            conn_comput_session_sub(:,:,:,k,kk) = source_conn{j}.(par);
            Seed_idx_session(k,kk) = datainfo.Seed_idx{j};
            k=k+1;
        end
    end
    conn_comput_session_sub_all{k1} = conn_comput_session_sub;
    Seed_idx_session_all{k1} = unique(Seed_idx_session);
    k1=k1+1;
end
size(conn_comput_session_sub_all{2})
Seed_idx_session_all{1}

% foi = [1,12]; foi = [8,12]; foi = [17,25]; foi = [2,7];
foi = [1,30];
for subsel = 1:length(datainfo.UQSub)
    for kk=1:size(conn_comput_session_sub_all{subsel},5)
        aedge = conn_comput_session_sub_all{subsel}(:,:,:,:,kk);
        size(aedge)
        aedge = mean(aedge,4);
        [~, idx_mn] = min(abs(foi(1) - source_conn{1}.freq));
        [~, idx_mx] = min(abs(foi(2) - source_conn{1}.freq));
        aedge = abs(aedge);
        aedge =  mean(aedge(:,:,idx_mn:idx_mx),3);
        seed_idx = Seed_idx_session_all{subsel};
        aedge(isnan(aedge))=0; tedge = (aedge.* double(aedge > 0.5.*max(aedge(:))));
        saveid = ['32037_00', num2str(sub_id(subsel)), ' session: ', num2str(kk)];
        close all
%         conn_mesh_display_vy('', '', '', Coordinate, tedge, 0.2);
        savepath = '/data/MEG/Research/logopenicppa/results/1_30/Conn_maps_wholebrain';
        savefig = fullfile(savepath, [saveid,'.jpg']);
        save(fullfile(savepath, [saveid,'_', conn_method, '.mat']), 'tedge');
%         conn_print(savefig,...
%             '-nogui',...
%             '-row',...
%             get(findobj(gcf,'label','Left view (both hem)'),'callback'),...
%             get(findobj(gcf,'label','Anterior view'),'callback'),...
%             get(findobj(gcf,'label','Superior view'),'callback'));
        pause (0.1),
        close all
    end
end

%- Progression analysis
d = rdir(fullfile(savepath, [conn_method, '*.mat']));

mm = [];
conn_mat_all = [];
for i=1:length(d)
    conn_mat = load(d(i).name);
    disp(d(i).name)
    conn_mat_all(i,:,:) = conn_mat.tedge;
    mm(i) = nanmean(nanmean(conn_mat.tedge(:)));
end

figure, bar(mm)

%%
% % d  = rdir(fullfile(outdir,'/**/Rest/**/seed*1-30Hz.mat'));
% for i=1:length(d)
%     %     conn_comput = load(d(i).name);
%     foi = conn_comput.net_conn_seed.foi;
%     source_conn = conn_comput.net_conn_seed.source_conn;
%     par = 'wpli_debiasedspctrm';
%     seed_idx = conn_comput.net_conn_seed.seed_idx;
%     
%     [~, idx_mn] = min(abs(foi(1) - source_conn.freq));
%     [~, idx_mx] = min(abs(foi(2) - source_conn.freq));
%     aedge =  mean(source_conn.(par)(:,:,idx_mn:idx_mx),3);
%     
%     tkz = tokenize(d(i).name,'/');
%     Sub{i} = tkz{end-3}; Session{i} = tkz{end-1}; run{i} = tkz{end}(end-11);
%     
%     aedge1 = zeros(size(aedge));
%     aedge1(seed_idx,:) = aedge(seed_idx,:);
%     aedge1(:,seed_idx) = aedge(:,seed_idx);
%     close all
%     conn_mesh_display_vy('', '', '', Coordinate, aedge1, 0.2);
%     %     %     pause,
%     disp(d(i).name)
%     disp([Sub{i}, '_', Session{i}, '_', run{i},'.jpg'])
%     %     pause(0.1),
%     savepath = '/data/MEG/Research/logopenicppa/results/1_30/Conn_maps';
%     savefig = fullfile(savepath, [Sub{i}, '_', Session{i}, '_', run{i},'.jpg']);
%     conn_print(savefig,...
%         '-nogui',...
%         '-row',...
%         get(findobj(gcf,'label','Left view (both hem)'),'callback'),...
%         get(findobj(gcf,'label','Anterior view'),'callback'),...
%         get(findobj(gcf,'label','Superior view'),'callback'));
%     pause(0.1),
%     close all
% end

%%
% % tmplate = load('temp_grid') % low-res
% % addpath('/data/MEG/Vahab/Github/MCW_MEGlab/FT/functions/External/brewermap')
% % 
% % 
% % seed_conn_sample = load('/data/MEG/Research/logopenicppa/ft_process/32037_002/Rest/1/seed_conn_32037_002_run_2_1-30Hz');
% % coor = seed_conn_sample.net_conn_seed.vs_roi.coor;
% 
% foi = conn.net_conn_seed.foi;
% source_conn = conn.net_conn_seed.source_conn;
% par = 'wpli_debiasedspctrm';
% seed_idx = conn.net_conn_seed.seed_idx;
% 
% n = 5;
% tedge = squeeze(conn.net_conn_seed.source_conn.wpli_debiasedspctrm(:,:,1));
% 
% [~, idx_mn] = min(abs(foi(1) - source_conn.freq));
% [~, idx_mx] = min(abs(foi(2) - source_conn.freq));
% % aedge =  squeeze(source_conn1.(par)(:,:,1));
% aedge =  mean(source_conn.(par)(:,:,idx_mn:idx_mx),3);
% % aedge(nonIdenticalIndices,:) = 0;
% % size(aedge)
% aedge1 = zeros(size(aedge));
% aedge1(seed_idx,:) = aedge(seed_idx,:);
% aedge1(:,seed_idx) = aedge(:,seed_idx);
% 
% 
% 
% conn_map = 1e3.*aedge(seed_idx,seed_idx);
% mni = (Coordinate(seed_idx,:));
% 
% % tedge= 1*rand(3,3);
% % mni = [1,20,10; -10,1,-20; 60,2,25];
% conn_mesh_display('', '', '', Coordinate, aedge1, .2);
% 
% figure,ft_plot_mesh(coor);

%%
template_grid.coordsys = 'mni';
cfg = [];
cfg.atlas      = atlas;
cfg.roi        = atlas.tissuelabel;  % here you can also specify a single label, i.e. single ROI
cfg.inputcoord = 'mni';
mask           = ft_volumelookup(cfg, template_grid);

%%
% Label = zeros(length(individual_grid.pos),1);

for jj=1:length(atlas.tissuelabel)
    
    cfg = [];
    cfg.atlas      = atlas;
    cfg.roi        = atlas.tissuelabel{jj};  % here you can also specify a single label, i.e. single ROI
    cfg.inputcoord = 'mni';
    mask1           = ft_volumelookup(cfg, template_grid);
    
    individual_grid3 = template_grid;
    individual_grid3.inside = false(individual_grid3.dim);
    individual_grid3.inside(mask1==1) = true;
    
    %                 hold on;
    %                 ft_plot_mesh(individual_grid3.pos(individual_grid3.inside,:),'vertexcolor', [rand(1,1),rand(1,1),rand(1,1)]);
    %             'edgecolor', [rand(1,1),rand(1,1),rand(1,1)]);
    grid2 = template_grid;
    grid2.inside = reshape(grid2.inside,length(grid2.pos),1);
    
    find(grid2.inside==1)
    Label(grid2.inside==1,:)= jj;
    
end
% figure, bar(Label(Label~=0))

%%
Coordinate = [
-48, -4, 44
50, -6, 44
-26, 56, 10
28, 56, 10
-16, 62, -12
18, 62, -12
-22, 10, 46
26, 10, 46
-12, 50, -8
12, 50, -8
-44, 24, 18
44, 24, 18
-34, 34, -2
34, 34, -2
-20, 52, -24
20, 52, -24
-50, -18, 20
52, -16, 22
-10, 14, 62
12, 14, 62
-12, 16, -18
10, 16, -18
-4, 52, 30
6, 52, 30
-12, 48, -10
14, 48, -10
-6, 38, -14
6, 38, -14
-38, 16, 0
38, 16, 0
-6, 24, 28
8, 24, 28
-2, 12, 26
4, 12, 26
-6, -4, 34
6, -4, 34
-26, -30, -8
28, -30, -8
-20, -22, -20
22, -22, -20
-18, -8, -20
20, -8, -20
-10, -74, 6
10, -74, 6
-12, -82, 6
12, -82, 6
-16, -64, -6
16, -64, -6
-18, -86, 10
18, -86, 10
-22, -90, 20
22, -90, 20
-26, -92, -6
26, -92, -6
-38, -52, -16
38, -52, -16
-50, -24, 38
50, -24, 38
-14, -46, 66
14, -46, 66
-42, -52, 44
42, -52, 44
-46, -40, 40
46, -40, 40
-46, -62, 34
46, -62, 34
-8, -68, 60
8, -68, 60
-6, -24, 60
6, -24, 60
-10, 12, 14
10, 12, 14
-22, 2, 2
26, 4, -4
-22, -4, -4
22, -4, -4
-12, -18, 8
12, -18, 8
-52, -18, 8
52, -18, 8
-56, -24, -8
56, -24, -8
-42, 14, -30
42, 14, -30
-60, -38, 6
60, -38, 6
-38, 14, -38
38, 14, -38
-46, -68, -12
46, -68, -12
-30, -68, -20
30, -68, -20
-18, -56, -16
18, -56, -16
-16, -68, -30
16, -68, -30
-22, -62, -46
22, -62, -46
-26, -44, -56
26, -44, -56
-34, -38, -60
34, -38, -60
-42, -30, -54
42, -30, -54
-50, -20, -44
50, -20, -44
-58, -12, -34
58, -12, -34
0, -64, -26
0, -56, -34
0, -48, -42
0, -38, -48
0, -30, -56
0, -20, -62
0, -12, -64
0, -4, -64];







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