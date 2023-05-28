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
addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/Logopenicppa/functions')
[atlas, all_path] = do_setpath();

%%
indir = '/data/MEG/Research/logopenicppa';
outdir = fullfile(indir,'ft_process');
rawdatadir = fullfile(indir,'raw');

%%
aal_coordinates

%%
conn_method = 'amplcorr'; 
conn_method = 'wpli_debiased';
par = [conn_method, 'spctrm'];

idx = [12:2:20,80:2:88];
roiidx = [idx, idx-1];

idx = [12:2:20,80:2:88];
idx(5:6)=[];
roiidx = [idx, idx-1];

label = {'Frontal-Inf-Oper-R' ;
    'Frontal-Inf-Tri-R'  ;
    'Frontal-Inf-Orb-R'  ;
    'Rolandic-Oper-R'    ;
    'Supp-Motor-Area-R'  ;
    'Heschl-R'           ;
    'Temporal-Sup-R'     ;
    'Temporal-Pole-Sup-R';
    'Temporal-Mid-R'     ;
    'Temporal-Pole-Mid-R';
    'Frontal-Inf-Oper-L' ;
    'Frontal-Inf-Tri-L'  ;
    'Frontal-Inf-Orb-L'  ;
    'Rolandic-Oper-L'    ;
    'Supp-Motor-Area-L'  ;
    'Heschl-L'           ;
    'Temporal-Sup-L'     ;
    'Temporal-Pole-Sup-L';
    'Temporal-Mid-L'     ;
    'Temporal-Pole-Mid-L'};

%%
clc
ft_progress('init', 'text',     'please wait ...');
d  = rdir(fullfile(outdir,['/**/Rest/**/sem*',conn_method, '.mat']));
clear Sub Run Session source_conn network network_atlas Seed_idx
for i=1:length(d)
    ft_progress(i/length(d), 'Processing network values %d from %d', i, length(d));
%     disp(d(i).name)
    tkz = tokenize(d(i).name,'/');
    conn_comput = load(d(i).name);
    network_atlas{i} = conn_comput.net_conn_seed.network_atlas;
    network{i} = conn_comput.net_conn_seed.network;
    source_conn{i} = conn_comput.net_conn_seed.source_conn;
    pause(0.1);
    Sub{i} = tkz{end-3}; Session{i} = tkz{end-1};
    tkz1 = tokenize(tkz{end},'_');
    run{i} = tkz1{end-1};
end
ft_progress('close');

%%
datainfo = [];
datainfo.Sub = Sub;
datainfo.Session = Session;
datainfo.run = run;
datainfo.UQSub = unique((datainfo.Sub));

%%
for i=1:1%length(d)
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
clc
close all
k=1; k1=1;
clear conn_sub_all
for i=[2,4,5,6]
    idx = find(strcmp(Sub, ['32037_00', num2str(i)])==1);
    disp(idx)
    k=1;
    clear conn_sub
    for j=idx
%         disp(j)
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

aedge =  mean(mm(:,:,idx_mn:idx_mx),3);

%% Map seed
% close all
% unq_seed_idx = unique(cell2mat(Seed_idx));
% for j=1:length(unq_seed_idx)
%     seed_map = atlas;
%     seed_map.tissue = zeros(size(atlas.tissue));
%     for i = 1:116
%         idx = atlas.tissue == unq_seed_idx(j);
%         seed_map.tissue(idx) = 1;
%     end
%     
%     cfg = [];
%     cfg.subj = ['seed:', atlas.tissuelabel{unq_seed_idx(j)}];
%     cfg.mask = 'tissue';
%     cfg.thre = 0;
%     cfg.savepath = '/data/MEG/Research/logopenicppa/results/1_30/Seed_roi';
%     cfg.colorbar = 0;
%     cfg.saveflag = 1;
%     cfg.colormap = [];
%     cfg.surfinflated   = 'surface_inflated_both.mat';
%     % cfg.views = [90 0;0,90;-90 0;-180,-90];
%     cfg.views = [-90,0; 90,0];
%     cfg.tit = ['seed:', atlas.tissuelabel{unq_seed_idx(j)}];
%     do_mapvis(cfg, seed_map);    
% end

%% Subj - group
% mnetwork_atlas = network_atlas{1};
close all
sub_id = [2,4,5,6];
k=1; k1=1;
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

for sel = 1:length(sub_id)
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
        cfg.savepath = '/data/MEG/Research/logopenicppa/results/1_30/Sem_conn_session';
        cfg.colorbar = 0;
        cfg.saveflag = 1;
        % cfg.colormap = colormap(flipud(pink));
        cfg.colormap = [];
        cfg.surfinflated   = 'surface_inflated_both.mat';
        %     cfg.views = [90 0;0,90;-90 0;-180,-90];
        cfg.views = [-90,0; 90,0];
        cfg.tit = ['32037_00', num2str(sub_id(sel)), ' session: ', num2str(kk)];
        do_mapvis(cfg, mnetwork_atlas);
        pause
%         (0.1),
        close all
    end
end

%% Conn toolbox
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
%             Seed_idx_session(k,kk) = datainfo.Seed_idx{j};
            k=k+1;
        end
    end
    conn_comput_session_sub_all{k1} = conn_comput_session_sub;
%     Seed_idx_session_all{k1} = unique(Seed_idx_session);
    k1=k1+1;
end
% size(conn_comput_session_sub_all{2})
% Seed_idx_session_all{1}

foi = [1,30];
% foi = [4,7];
% foi = [18,25];
% foi = [8,12];
% foi = [4,30];
for subsel = 1:length(datainfo.UQSub)
    for kk=1:size(conn_comput_session_sub_all{subsel},5)        
        aedge = conn_comput_session_sub_all{subsel}(:,:,:,:,kk);
        size(aedge)
        aedge = mean(aedge,4);
        [~, idx_mn] = min(abs(foi(1) - source_conn{1}.freq));
        [~, idx_mx] = min(abs(foi(2) - source_conn{1}.freq));
%         aedge = abs(aedge);
        aedge =  mean(aedge(:,:,idx_mn:idx_mx),3);
        %         seed_idx = Seed_idx_session_all{subsel};
        saveid = ['32037_00', num2str(sub_id(subsel)), ' session: ', num2str(kk)];

%         aedge1 = zeros(size(aedge));
%         aedge1(seed_idx,:) = aedge(seed_idx,:);
%         aedge1(:,seed_idx) = aedge(:,seed_idx);
        
        aedge1 = zeros(size(aedge));
        aedge1(roiidx,:) = aedge(roiidx,:);
        aedge1(:,roiidx) = aedge(:,roiidx);
        
        close all
        %         conn_mesh_display_vy('', '', '', Coordinate, aedge1, 0.2);
        savepath = '/data/MEG/Research/logopenicppa/results/1_30/Semantic_conn_maps';
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

% addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/FT_fucntions/External/other')
% mm = [];
% net_mat_all = [];
% for i=1:length(d)
%     conn_mat = load(d(i).name);
%     disp(d(i).name);
%     tmp = conn_mat.aedge1;
%     tmp(isnan(tmp(:)))=0;
%     net_mat_all(i,:,:) = eigenvector_centrality_und(tmp);
%     mm(i) = nanmean(net_mat_all(i,:,:));
% end
% figure, bar(mm)


%%
temp = squeeze(mean(conn_mat_all(1:4,:,:),1)); size(temp)
conn_diff = squeeze(conn_mat_all(2,:,:) - conn_mat_all(1,:,:)); 
conn_diff = squeeze(conn_mat_all(4,:,:) - conn_mat_all(3,:,:)); 
conn_diff = squeeze(conn_mat_all(4,:,:) - temp); size(conn_diff)

conn_diff = squeeze(conn_mat_all(1,:,:));
size(conn_diff)
conn_diff (abs(conn_diff) < 0.5.*max(abs(conn_diff(:)))) = 0;
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

