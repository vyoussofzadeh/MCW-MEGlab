function [net_sel_mutiple_label, LI_sub] = do_group_LI_net_contrast1(cfg_main)

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

%%
% net_sel_mutiple_label = {'Angular'; 'Frontal'; 'Occipital'; 'Other'; 'PCingPrecun';'Temporal'; 'BTLA'; 'VWFA'};
net_sel_mutiple_label = Data_hcp_atlas.groups_labels';
% net_sel_mutiple_label = {'Angular'; 'Frontal'; 'Occipital'; 'Other'; 'PCingPrecun';'Temporal'; 'Whole'};
% idx_L{end+1} = idx_L_whole;idx_R{end+1} = idx_R_whole;

% S_data_sel{1}.subcon
% S_data_sel{2}.subcon


cd(data_save_dir)
savefilename = fullfile(data_save_dir,['LI_', S_data_sel{1}.subcon, '.mat']);

if exist(savefilename,'file') == 2
    load(savefilename),
else
    
    ft_progress('init', 'text',     'please wait ...');
    clear m_LI_max_sub LI_sub
    for j=1:length(net_sel_mutiple_label)
        ft_progress(j/length(net_sel_mutiple_label), 'Processing networks %d from %d', j, length(net_sel_mutiple_label));
        
        %- Li Calc.
        cfg = [];
        cfg.sinput = S_data_sel; %sFiles_anim_hc;
        cfg.BS_data_dir = BS_data_dir;
        cfg.atlas = Data_hcp_atlas.atlas; cfg.thre = thre; cfg.fplot = 0;
        cfg.index_L = idx_L{j};
        cfg.index_R = idx_R{j};
        cfg.Threshtype = Threshtype;
        
        for i=1:length(sFiles_in{1}.sFiles_in)
            pause(0.1);
            
            s_in = [sFiles_in{1}.sFiles_in(i); sFiles_in{2}.sFiles_in(i)];
            cfg.sinput = s_in;
            cfg.wi = wi;
            %             [LI, wi_max, pow] = do_lat_analysis_asymetric(cfg);
            switch method
                case 'wDICS_bsl_threshold'
                    [LI, ~, pow] = do_lat_analysis_asymetric_contrast(cfg);
                    pow_sub(j,i,:) = pow;
                case { 'DICS_baseline_Counting', 'LCMV_basline_Counting'}
                    [LI, ~] = do_lat_analysis_DICS_baseline_Counting(cfg);
                case 'wDICS_bsl_bootstrapping'
                    cfg.divs = 150;
                    cfg.n_resampling = 200;
                    cfg.RESAMPLE_RATIO = 0.9;
                    cfg.dvd = 5; % 5
                    [LI, ~] = do_LI_bootstrap_contrast(cfg);
            end
            %             LI_count = LI;
            %             close all, figure, plot(LI_count), hold on, plot(LI, 'r')
%             figure, plot(LI_count), hold on, plot(LI, 'r')
            
            LI_sub(j,i,:) = LI;
            m_LI_max_sub(i) = nanmean(LI);
        end
        
    end
    ft_progress('close')
    switch method
        case 'threshold'
            save(savefilename,'LI_sub','m_LI_max_sub','pow_sub'),
        otherwise
            save(savefilename,'LI_sub','m_LI_max_sub'),
    end
end

end