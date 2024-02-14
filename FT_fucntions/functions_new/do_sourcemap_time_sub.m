function ssinput = do_sourcemap_time_sub(cfg_main)

S_data_sel = cfg_main.S_data_sel;
wi_sub_max = cfg_main.wi_sub_max;
BS_data_dir = cfg_main.BS_data_dir;
sFiles_in = cfg_main.S_data_sel.sFiles_in;
src = cfg_main.src;

%%
subj_sel = input(['sel subject (out of)', num2str(length(sFiles_in)),':']);
ssinput = load(fullfile(BS_data_dir,sFiles_in{subj_sel}));
grand_val = ssinput.ImageGridAmp;

% toi_pst = L_max_hc_anim;
disp(wi_sub_max(subj_sel,:))
toi_pst = input('enter time of interests, e.g.,[0,2]: ');
Time = ssinput.Time;

[~, idx_toi1] = min(abs(Time - toi_pst(1))); 
[~, idx_toi2] = min(abs(Time - toi_pst(2)));

tmp  = mean(grand_val(:,idx_toi1:idx_toi2),2);
% tmp(abs(tmp) < thre.*max(tmp(:))) = 0;
% tmp = (tmp - min(tmp(:))) ./ (max(tmp(:)) - min(tmp(:))); %
% tmp = tmp./max(tmp(:));
%     tmp = tmp.^0.25;

figure,
cfg = [];
cfg.view = [-180,-90;0,90;-90,0; 90,0];
cfg.color = (viridis(256));
cfg.position = [800   800   1000   300];
cfg.title = S_data_sel.sFiles_subid{subj_sel}; 
cfg.alpha = 1; cfg.coor = [];
cfg.surf = src; cfg.d_in = tmp;
do_surfplot(cfg);

end