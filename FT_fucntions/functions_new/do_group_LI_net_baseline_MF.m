function [net_sel_mutiple_label, LI_sub] = do_group_LI_net_baseline_MF(cfg_main)

% m_LI_max_sub, pow_sub

sFiles_in = cfg_main.S_data_sel;
BS_data_dir = cfg_main.BS_data_dir;
S_data_sel = cfg_main.S_data_sel;
Data_hcp_atlas = cfg_main.Data_hcp_atlas;
thre = cfg_main.thre;
idx_L = cfg_main.idx_L;
idx_R = cfg_main.idx_R;
wi = cfg_main.wi;
data_save_dir = cfg_main.data_save_dir;
method = cfg_main.method;
Threshtype = cfg_main.Threshtype;
doavg = cfg_main.doavg;

%%
% net_sel_mutiple_label = {'Angular'; 'Frontal'; 'Occipital'; 'Other'; 'PCingPrecun';'Temporal'; 'BTLA'; 'VWFA'};
net_sel_mutiple_label = Data_hcp_atlas.groups_labels';
% net_sel_mutiple_label = {'Angular'; 'Frontal'; 'Occipital'; 'Other'; 'PCingPrecun';'Temporal'; 'Whole'};
% idx_L{end+1} = idx_L_whole;idx_R{end+1} = idx_R_whole;

% S_data_sel{1}.subcon
% S_data_sel{2}.subcon

%% Create a folder
if ~exist(data_save_dir, 'dir')
    mkdir(data_save_dir);
    disp('Folder created successfully.');
else
    disp('Folder already exists.');
end
cd(data_save_dir)

%%
savefilename = fullfile(data_save_dir,['LI_', S_data_sel{1}.subcon, '.mat']);

sel_netrois = [1,2,6,11];
if exist(savefilename,'file') == 2
    load(savefilename),
else
    
    ft_progress('init', 'text',     'please wait ...');
    clear m_LI_max_sub LI_sub
    for j = sel_netrois %1:length(net_sel_mutiple_label)
        
        ft_progress(j/length(net_sel_mutiple_label), 'Processing networks %d from %d', j, length(net_sel_mutiple_label));
        
        %- Li Calc.
        cfg = [];
        cfg.sinput = S_data_sel; %sFiles_anim_hc;
        cfg.BS_data_dir = BS_data_dir;
        cfg.atlas = Data_hcp_atlas.atlas;
        cfg.thre = thre;
        cfg.fplot = 0; % cfg.fplot = 1;
        cfg.index_L = idx_L{j};
        cfg.index_R = idx_R{j};
        cfg.Threshtype = Threshtype;
        cfg.doavg = doavg;
        cfg.parcellaion = 1;
        cfg.math = cfg_main.math; % 'db', ,'diff', 'rec_diff', 'rec_diff_mbsl', 'rec_diff_zbsl'
        
        for i=1:length(sFiles_in{1}.sFiles_in)
            pause(0.1);
            
            s_in = [sFiles_in{1}.sFiles_in(i); sFiles_in{2}.sFiles_in(i)];
            cfg.sinput = s_in;
            cfg.wi = wi;
            if contains(method, 'Magnitude')                
                [LI, ~, pow] = do_lat_analysis_asymetric_magnitude2_MF(cfg);
                pow_sub(j,i,:) = pow;
            elseif contains(method, 'Counting')
                [LI, ~, roi_count] = do_lat_analysis_baseline_Counting_MF(cfg);
                count_sub(j,i,:) = roi_count;
            elseif contains (method, 'Bootstrapping')
                cfg.divs = 25;
                cfg.n_resampling = 50;
                cfg.RESAMPLE_RATIO = 0.75;
                cfg.dvd = 5; % 5
                cfg.downsamplerate = 5; % 2 times down-sampling - has to be fixed!
                [LI, ~, ~, roi_count] = do_LI_bootstrap_contrast_MF(cfg);
                count_sub(j,i,:) = roi_count;
            end
            
            LI_sub(j,i,:) = LI;
            m_LI_max_sub(i) = nanmean(LI);
        end
        
    end
    ft_progress('close')
    
    setup = [];
    setup.BS_data_dir = BS_data_dir;
    setup.thre = thre;
    setup.Threshtype = Threshtype;
    setup.S_data = S_data_sel;
    setup.math = cfg.math;
    
    switch contains(method,'Magnitude')
        case 1
            save(savefilename,'LI_sub','m_LI_max_sub','pow_sub', 'setup'),
        otherwise
            save(savefilename,'LI_sub','m_LI_max_sub','count_sub', 'setup'),
    end
end
end