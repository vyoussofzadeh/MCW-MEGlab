function [groups_labels, groups] = do_read_HCP_labels_bs()


% link: https://github.com/mne-tools/mne-python/blob/main/mne/datasets/utils.py


groups = {{'Primary Visual Cortex (V1)','V1'}, ...
    {'Early Visual Cortex','V2', 'V3', 'V4'}, ...
    {'Dorsal Stream Visual Cortex','V3A', 'V3B', 'V6', 'V6A', 'V7', 'IPS1'}, ...
    {'Ventral Stream Visual Cortex','V8', 'VVC', 'PIT', 'FFC', 'VMV1', 'VMV2', 'VMV3'}, ...
    {'MT+ Complex and Neighboring Visual Areas', 'V3CD', 'LO1', 'LO2', 'LO3', 'V4t', 'FST', 'MT', 'MST', 'PH'},...
    {'Somatosensory and Motor Cortex','4', '3a', '3b', '1', '2'},...
    {'Paracentral Lobular and Mid Cingulate Cortex','24dd', '24dv', '6mp', '6ma', 'SCEF', '5m', '5L', '5mv'},...
    {'Premotor Cortex','55b', '6d', '6a', 'FEF', '6v', '6r', 'PEF'},...
    {'Posterior Opercular Cortex','43', 'FOP1', 'OP4', 'OP1', 'OP2-3', 'PFcm'},...
    {'Early Auditory Cortex','A1', 'LBelt', 'MBelt', 'PBelt', 'RI'},...
    {'Auditory Association Cortex','A4', 'A5', 'STSdp', 'STSda', 'STSvp', 'STSva', 'STGa', 'TA2'},...
    {'Insular and Frontal Opercular Cortex','52', 'PI', 'Ig', 'PoI1', 'PoI2', 'FOP2', 'FOP3', ...
    'MI', 'AVI', 'AAIC', 'Pir', 'FOP4', 'FOP5'},...
    {'Medial Temporal Cortex','H', 'PreS', 'EC', 'PeEc', 'PHA1', 'PHA2', 'PHA3'},...
    {'Lateral Temporal Cortex','PHT', 'TE1p', 'TE1m', 'TE1a', 'TE2p', 'TE2a', ...
    'TGv', 'TGd', 'TF'},...
    {'Temporo-Parieto-Occipital Junction','TPOJ1', 'TPOJ2', 'TPOJ3', 'STV', 'PSL'},...
    {'Superior Parietal Cortex','LIPv', 'LIPd', 'VIP', 'AIP', 'MIP', ...
    '7PC', '7AL', '7Am', '7PL', '7Pm'},...
    {'Inferior Parietal Cortex','PGp', 'PGs', 'PGi', 'PFm', 'PF', 'PFt', 'PFop', ...
    'IP0', 'IP1', 'IP2'},...
    {'Posterior Cingulate Cortex','DVT', 'ProS', 'POS1', 'POS2', 'RSC', 'v23ab', 'd23ab', ...
    '31pv', '31pd', '31a', '23d', '23c', 'PCV', '7m'},...
    {'Anterior Cingulate and Medial Prefrontal Cortex','33pr', 'p24pr', 'a24pr', 'p24', 'a24', 'p32pr', 'a32pr', 'd32', ...
    'p32', 's32', '8BM', '9m', '10v', '10r', '25',},...
    {'Orbital and Polar Frontal Cortex','47s', '47m', 'a47r', '11l', '13l',...
    'a10p', 'p10p', '10pp', '10d', 'OFC', 'pOFC',},...
    {'Inferior Frontal Cortex','44', '45', 'IFJp', 'IFJa', 'IFSp', 'IFSa', '47l', 'p47r',},...
    {'DorsoLateral Prefrontal Cortex','8C', '8Av', 'i6-8', 's6-8', 'SFL', '8BL', '9p', '9a', '8Ad',...
    'p9-46v', 'a9-46v', '46', '9-46d',}};


%%
groups_labels = {{'Prim Visual C. (V1)'}, ...
    {'Early Vis C.'}, ...
    {'Dors Stream Vis C.'}, ...
    {'Vent Stream Vis C.'}, ...
    {'MT+ Complex and Neighboring Vis Areas'},...
    {'Somatosensory and Motor C.'},...
    {'Paracentral Lobular and Mid Cin C.',},...
    {'Premotor C.'},...
    {'Post Opercular C.'},...
    {'Early Aud C.'},...
    {'Aud Association C.'},...
    {'Insular and Front Oper C.'},...
    {'Med Temp C.'},...
    {'Lateral Temp C.'},...
    {'Temporo-Parieto-Occip Junction'},...
    {'Sup Pari C.'},...
    {'Inf Pari C.'},...
    {'Post Cin C.'},...
    {'Ant Cin and Med PreFront C.'},...
    {'Orbital and Polar Front C.'},...
    {'Inf Front C.'},...
    {'DorsoLateral PreFront C.'}};


end
%%
% idx_L32 = [];
% for i=1:length(groups)
%     grois = groups{i}(2:end);
%     idx = [];
%     for j=1:length(grois)
%             idx(j) = strmatch(['L_',grois{j}, '_ROI'],cfg.rois);
%     end
%     idx_L32{i} = idx;
% end
% 
% idx_R32 = [];
% for i=1:length(groups)
%     grois = groups{i}(2:end);
%     idx = [];
%     for j=1:length(grois)
%             idx(j) = strmatch(['R_',grois{j}, '_ROI'],cfg.rois);
%     end
%     idx_R32{i} = idx;
% end
% 
% 
% %%
% % idx_L = cfg.lat_index(:,1);
% % idx_R = cfg.lat_index(:,2);
% % idx_lr = [idx_L;idx_R];
% % 
% % cfg = [];
% % cfg.atlas = atlas;
% % cfg.src_fname = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/data_info/cortex_pial_low.fs';
% % cfg.sel = 'roi'; % 'whole', 'left', 'right', 'roi';
% % cfg.lat_index = [idx_L, idx_R];
% % cfg.rois = rois;
% % do_plot_atlas(cfg)
% 
% % % close all
% % idx_LR32 = [idx_L32,idx_R32];
% % length(idx_lr32)
% % length(idx_L)
% % length(idx_R)
% 
% %- Whole atlas
% vertexcolor = zeros(size(src.pos,1), 3);
% for iScout=1:nScouts
%     index = Scouts(iScout).Vertices;
%     if ~isempty(index)
%         vertexcolor(index,:) = repmat(Scouts(iScout).Color,  length(index), 1);
%         vertexcolor(index,:) = repmat([0.5,0.5,0.5],  length(index), 1);
%     end
% end
% 
% 
% % left ROIs
% % for iScout=1:length(idx_L32)
% %     for j=1:length(idx_L32{iScout})
% %         index = Scouts((idx_L32{iScout}(j))).Vertices;
% %         if ~isempty(index)
% %             vertexcolor(index,:) = repmat(Scouts(idx_lr(iScout)).Color,  length(index), 1);
% %         end
% %     end
% % end
% 
% % % right ROIs
% % for iScout=1:length(idx_R32)
% %     for j=1:length(idx_R32{iScout})
% %         index = Scouts((idx_R32{iScout}(j))).Vertices;
% %         if ~isempty(index)
% %             vertexcolor(index,:) = repmat(Scouts(idx_lr(iScout)).Color,  length(index), 1);
% %         end
% %     end
% % end
% 
% 
% % close all
% % left_sel ROIs
% groups_labels_num = [];
% for i=1:length(groups_labels)
%     groups_labels_num{i} = [num2str(i), ': ', groups_labels{i}{1}];
% end
% disp(cell2table(groups_labels_num'));
% sel = input('enter rois:');
% for iScout=1:length(sel)
%     for j=1:length(idx_L32{sel(iScout)})
%         index = Scouts((idx_L32{sel(iScout)}(j))).Vertices;
%         if ~isempty(index)
%             vertexcolor(index,:) = repmat(Scouts(idx_lr(iScout)).Color,  length(index), 1);
%         end
%     end
% end
% 
% 
% % Visualisation de l'atlas Desikan_killiany
% % close all
% % figure;
% % ft_plot_mesh(src, 'faecolor', 'brain',  'vertexcolor', ...
% %     vertexcolor, 'facealpha', 1);
% % view(-3, 2);
% 
% % addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/FT_fucntions/External/Colormaps-from-MatPlotLib2.0')
% 
% figure
% cfg = [];
% cfg.view = [-180,-90;0,90;-90,0; 90,0; 0, 0];
% cfg.position = [800   800   1000   300];
% cfg.color = (viridis(256));
% cfg.title = ['']; cfg.alpha = 1; cfg.coor = [];
% do_surfplot(cfg,src,vertexcolor);
% % title([num2str(sel), ': ', rois{idx_L_Sel}])


%%



% end