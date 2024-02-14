function m_source = do_sourcemap_time_sub_optimal_toi_diff(cfg_main)

% clc
% src_fname = cfg_main.src_fname;
% Data_hcp_atlas = cfg_main.Data_hcp_atlas;
% network_sel = cfg_main.network_sel;

% close all
% figure,
% ssinput = cfg_main.ssinput;

sFiles_in1 = cfg_main.S_data_sel1.sFiles_in;
sFiles_in2 = cfg_main.S_data_sel2.sFiles_in;

BS_data_dir = cfg_main.BS_data_dir;
wi_sub_max = cfg_main.wi_sub_max;
src = cfg_main.src;

clear s_all
ft_progress('init', 'text',     'please wait ...');
for subj_sel=1:length(sFiles_in1)
    
    ssinput1 = load(fullfile(BS_data_dir,sFiles_in1{subj_sel}));
    ssinput2 = load(fullfile(BS_data_dir,sFiles_in2{subj_sel}));
    
    Time = ssinput1.Time;
    ft_progress(subj_sel/length(sFiles_in1), 'Processing subjects %d from %d', subj_sel, length(sFiles_in1));
    grand_val = ssinput1.ImageGridAmp - ssinput2.ImageGridAmp;
    toi_pst = wi_sub_max(subj_sel,:);
    
    [~, idx_toi1] = min(abs(Time - toi_pst(1)));
    [~, idx_toi2] = min(abs(Time - toi_pst(2)));
    
    tmp  = mean(grand_val(:,idx_toi1:idx_toi2),2);
    tmp = (tmp - min(tmp(:))) ./ (max(tmp(:)) - min(tmp(:)));
    
    %     cfg = [];
    %     cfg.view = [-180,-90;0,90;-90,0; 90,0];
    %     cfg.color = (viridis(256));
    %     cfg.position = [800   800   1000   300];
    %     cfg.title = sub_sel{subj_sel};
    %     cfg.alpha = 1; cfg.coor = [];
    %     cfg.surf = src; cfg.d_in = tmp;
    %     do_surfplot(cfg);
    %     pause(0.2)
    
    s_all(subj_sel,:) =  tmp;
end
ft_progress('close')

figure,
cfg = [];
cfg.view = [-180,-90;0,90;-90,0; 90,0];
cfg.color = viridis(256);
cfg.position = [800   800   1000   300];
cfg.title = '';
cfg.alpha = 1; cfg.coor = [];
m_source = mean(s_all,1)';
% tmp = (tmp - min(tmp(:))) ./ (max(tmp(:)) - min(tmp(:)));
cfg.surf = src; cfg.d_in = m_source;
do_surfplot(cfg);

end