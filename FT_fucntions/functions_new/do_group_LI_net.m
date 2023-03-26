function [net_sel_mutiple_label, LI_sub, m_LI_max_sub] = do_group_LI_net(cfg_main)

sFiles_in = cfg_main.S_data_sel.sFiles_in;
BS_data_dir = cfg_main.BS_data_dir;
S_data_sel = cfg_main.S_data_sel; 
Data_hcp_atlas = cfg_main.Data_hcp_atlas;
thre = cfg_main.thre;
idx_L = cfg_main.idx_L;
idx_R = cfg_main.idx_R;
wi = cfg_main.wi;
data_save_dir = cfg_main.data_save_dir;

%%
net_sel_mutiple_label = {'Angular'; 'Frontal'; 'Occipital'; 'Other'; 'PCingPrecun';'Temporal'; 'BTLA'; 'VWFA'};
% net_sel_mutiple_label = {'Angular'; 'Frontal'; 'Occipital'; 'Other'; 'PCingPrecun';'Temporal'; 'Whole'};
% idx_L{end+1} = idx_L_whole;idx_R{end+1} = idx_R_whole;

cd(data_save_dir)
savefilename = fullfile(data_save_dir,['LI_',S_data_sel.s_tag, '.mat']);
if exist(savefilename,'file') == 2
    load(savefilename),
else
    
    ft_progress('init', 'text',     'please wait ...');
    clear m_LI_max_sub LI_sub
    for j=1:length(net_sel_mutiple_label)
        ft_progress(j/length(net_sel_mutiple_label), 'Processing networks %d from %d', j, length(net_sel_mutiple_label));
        
        %- Li Calc.
        cfg = [];
        cfg.sinput = S_data_sel.sFiles_in; %sFiles_anim_hc;
        cfg.BS_data_dir = BS_data_dir;
        cfg.atlas = Data_hcp_atlas.atlas; cfg.thre = thre; cfg.fplot = 0;
        cfg.index_L = idx_L{j};
        cfg.index_R = idx_R{j};
        
        
        for i=1:length(sFiles_in)
            pause(0.1);
            cfg.sinput = sFiles_in{i};
            %     cfg.wi = wi_sub_max(i,:);
            cfg.wi = wi;
            [LI, wi_max] = do_lat_analysis_asymetric(cfg);
            LI_sub(j,i,:) = LI;
            m_LI_max_sub(i) = nanmean(LI);
            %     wi_sub_max(i,:) = wi_max;
        end
        
    end
    ft_progress('close')
%     save(savefilename,'Comment','LI_sub','m_LI_max_sub'),
    save(savefilename,'LI_sub','m_LI_max_sub'),
end

end