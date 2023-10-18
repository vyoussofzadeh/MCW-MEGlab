function Data_hcp_atlas = ecpfunc_hcp_atlas2(cfg_main)
% ECP functions
% Project: ECP_SD
% Written by: Vahab Youssof Zadeh
% Update: 05/31/2023

src_fname = cfg_main.src_fname;
glass_dir = cfg_main.glass_dir;
glass_atlas = cfg_main.glass_atlas;

%%

% '/group/jbinder/ECP/MEG/laterality_index/bilateral_glasser_lateral.tsv'

%%
% Load atlas
% atlas = load('/data/MEG/Vahab/Github/MCW_MEGlab/tools/Atlas/HCP/HCP atlas for Brainstorm/Best/scout_mmp_in_mni_symmetrical_final_updated.mat');
atlas = load(glass_atlas);
groups_labels = {'Angular', 'Frontal', 'Occipital', 'Other', 'PCingPrecun', 'Temporal'};

rois = {atlas.Scouts.Label};
region = {atlas.Scouts.Region};

all_idx_L = find(startsWith(rois, 'L_'))';
all_idx_R = find(startsWith(rois, 'R_'))';

if cfg_main.plotflag == 1;
    
    % Plot HCP atlas
    cfg = struct();
    cfg.atlas = atlas;
    cfg.src_fname = src_fname;
    % '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/Atlas/cortex_pial_low.fs';
    cfg.sel = 'roi'; % 'whole', 'left', 'right', 'roi';
    cfg.index_L = all_idx_L;
    cfg.index_R = all_idx_R;
    cfg.rois = rois;
    cfg.rois_sel = 1:180;
    cfg.title = '';
    do_plot_HCP_atlas(cfg)
    
end

% Update frontal region
idx_sel_L = strcmp(region(all_idx_L), 'LF');
idx_sel_R = strcmp(region(all_idx_R), 'RF');

% Load atlas ROIs from fMRI study
% '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/Atlas/Glasser';
load(fullfile(glass_dir, 'LI_glasser_manual_net_12.mat'), 'glass_net_L_label', 'glass_net_R_label');

% Update frontal region labels
glass_net_L_label{2} = [glass_net_L_label{2}; rois(all_idx_L(idx_sel_L))'];
glass_net_R_label{2} = [glass_net_R_label{2}; rois(all_idx_R(idx_sel_R))'];

% Add BTLA labels
btla = [2, 3, 5, 8, 9, 16, 17, 18, 21, 22]; net_sel = 6;
% glass_net_L_label{6} = glass_net_L_label{6}(btla);
% glass_net_R_label{6} = glass_net_R_label{6}(btla);

BTLA_L_label = [];
BTLA_R_label = [];
for i=1:length(btla)
    BTLA_L_label = [BTLA_L_label; glass_net_L_label{net_sel}(btla(i))];
    BTLA_R_label = [BTLA_R_label; glass_net_R_label{net_sel}(btla(i))];
end

glass_net_L_label{7} = BTLA_L_label;
glass_net_R_label{7} = BTLA_R_label;

groups_labels{7} = 'BTLA';

%% Add VWFA labels
vw2 = [6, 14, 15, 81]; net_sel = 4;
% glass_net_L_label{4} = glass_net_L_label{4}(vw2);
% glass_net_R_label{4} = glass_net_R_label{4}(vw2);

VW_L_label = [];
VW_R_label = [];
for i=1:length(vw2)
    VW_L_label = [VW_L_label; glass_net_L_label{net_sel}(vw2(i))];
    VW_R_label = [VW_R_label; glass_net_R_label{net_sel}(vw2(i))];
end

glass_net_L_label{8} = VW_L_label;
glass_net_R_label{8} = VW_R_label;

groups_labels{8} = 'VWFA';

% Add LT and RT region labels
% idx_sel_L = strcmp(region(all_idx_L), 'LT');
% idx_sel_R = strcmp(region(all_idx_R), 'RT');
% glass_net_L_label{7} = rois(all_idx_L(idx_sel_L));
% glass_net_R_label{7} = rois(all_idx_R(idx_sel_R));


%% Add ATG labels
ATG_labels = {'L_TGv_ROI', 'L_TGd_ROI'};

% Find indices that correspond to ATG and STG labels
ATG_L_label = [];
ATG_R_label = [];
for i=1:length(ATG_labels)
    ATG_L_label = [ATG_L_label; rois(strcmp(rois, ATG_labels{i}))];
    ATG_R_label = [ATG_R_label; rois(strcmp(rois, strrep(ATG_labels{i}, 'L_', 'R_')))];
end

glass_net_L_label{9} = ATG_L_label;
glass_net_R_label{9} = ATG_R_label;

groups_labels{9} = 'ATG'; %Anterior temporal G.

%% Post. STG
PSTG_labels = {'L_TE1p_ROI'};

% Find indices that correspond to ATG and STG labels
PSTG_L_label = [];
PSTG_R_label = [];
for i=1:length(PSTG_labels)
    PSTG_L_label = [PSTG_L_label; rois(strcmp(rois, PSTG_labels{i}))];
    PSTG_R_label = [PSTG_R_label; rois(strcmp(rois, strrep(PSTG_labels{i}, 'L_', 'R_')))];
end

glass_net_L_label{10} = PSTG_L_label;
glass_net_R_label{10} = PSTG_R_label;

groups_labels{10} = 'PSTG'; %Post. STG

%%
% ATG_L_label = rois(all_idx_L(ATG_indices));
% ATG_R_label = rois(all_idx_R(ATG_indices));
%
% glass_net_L_label{9} = ATG_L_label;
% glass_net_R_label{9} = ATG_R_label;
% groups_labels{9} = 'ATG'; %Anterior temporal G.
%
% % Add STG labels
% STG_indices = ... % Fill this with indices or identifiers for the STG region
% STG_L_label = rois(all_idx_L(STG_indices));
% STG_R_label = rois(all_idx_R(STG_indices));
%
% glass_net_L_label{10} = STG_L_label;
% glass_net_R_label{10} = STG_R_label;
% groups_labels{10} = 'STG'; %Superior temporal G.

%% Lateral rois
% see, /data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/Pipe_check_atlas.m for details
LI_glasser_lateral_rois = load(fullfile(glass_dir,'LI_glasser_lateral_rois.mat'));
groups_labels = [groups_labels, 'lateral'];

glass_net_L_label{11} = LI_glasser_lateral_rois.glass_roi_lat_L_name;
glass_net_R_label{11} = LI_glasser_lateral_rois.glass_roi_lat_R_name;


%%
Data_hcp_atlas.glass_net_L_label = glass_net_L_label;
Data_hcp_atlas.glass_net_R_label = glass_net_R_label;
Data_hcp_atlas.groups_labels = groups_labels;
Data_hcp_atlas.atlas = atlas;
Data_hcp_atlas.rois = rois;

end

% roi =     {'L_V1_ROI'    }
%     {'R_V1_ROI'    }
%     {'L_FEF_ROI'   }
%     {'R_FEF_ROI'   }
%     {'L_OP4_ROI'   }
%     {'R_OP4_ROI'   }
%     {'L_OP1_ROI'   }
%     {'R_OP1_ROI'   }
%     {'L_OP2-3_ROI' }
%     {'R_OP2-3_ROI' }
%     {'L_52_ROI'    }
%     {'R_52_ROI'    }
%     {'L_RI_ROI'    }
%     {'R_RI_ROI'    }
%     {'L_PFcm_ROI'  }
%     {'R_PFcm_ROI'  }
%     {'L_PoI2_ROI'  }
%     {'R_PoI2_ROI'  }
%     {'L_TA2_ROI'   }
%     {'R_TA2_ROI'   }
%     {'L_FOP4_ROI'  }
%     {'R_FOP4_ROI'  }
%     {'L_MI_ROI'    }
%     {'R_MI_ROI'    }
%     {'L_PEF_ROI'   }
%     {'R_PEF_ROI'   }
%     {'L_Pir_ROI'   }
%     {'R_Pir_ROI'   }
%     {'L_AVI_ROI'   }
%     {'R_AVI_ROI'   }
%     {'L_AAIC_ROI'  }
%     {'R_AAIC_ROI'  }
%     {'L_FOP1_ROI'  }
%     {'R_FOP1_ROI'  }
%     {'L_FOP3_ROI'  }
%     {'R_FOP3_ROI'  }
%     {'L_FOP2_ROI'  }
%     {'R_FOP2_ROI'  }
%     {'L_PFt_ROI'   }
%     {'R_PFt_ROI'   }
%     {'L_AIP_ROI'   }
%     {'R_AIP_ROI'   }
%     {'L_EC_ROI'    }
%     {'R_EC_ROI'    }
%     {'L_PreS_ROI'  }
%     {'R_PreS_ROI'  }
%     {'L_55b_ROI'   }
%     {'R_55b_ROI'   }
%     {'L_H_ROI'     }
%     {'R_H_ROI'     }
%     {'L_ProS_ROI'  }
%     {'R_ProS_ROI'  }
%     {'L_PeEc_ROI'  }
%     {'R_PeEc_ROI'  }
%     {'L_STGa_ROI'  }
%     {'R_STGa_ROI'  }
%     {'L_PBelt_ROI' }
%     {'R_PBelt_ROI' }
%     {'L_A5_ROI'    }
%     {'R_A5_ROI'    }
%     {'L_PHA1_ROI'  }
%     {'R_PHA1_ROI'  }
%     {'L_PHA3_ROI'  }
%     {'R_PHA3_ROI'  }
%     {'L_STSda_ROI' }
%     {'R_STSda_ROI' }
%     {'L_STSdp_ROI' }
%     {'R_STSdp_ROI' }
%     {'L_V3A_ROI'   }
%     {'R_V3A_ROI'   }
%     {'L_STSvp_ROI' }
%     {'R_STSvp_ROI' }
%     {'L_TGd_ROI'   }
%     {'R_TGd_ROI'   }
%     {'L_TE1a_ROI'  }
%     {'R_TE1a_ROI'  }
%     {'L_TE1p_ROI'  }
%     {'R_TE1p_ROI'  }
%     {'L_TE2a_ROI'  }
%     {'R_TE2a_ROI'  }
%     {'L_TF_ROI'    }
%     {'R_TF_ROI'    }
%     {'L_TE2p_ROI'  }
%     {'R_TE2p_ROI'  }
%     {'L_PHT_ROI'   }
%     {'R_PHT_ROI'   }
%     {'L_PH_ROI'    }
%     {'R_PH_ROI'    }
%     {'L_TPOJ1_ROI' }
%     {'R_TPOJ1_ROI' }
%     {'L_RSC_ROI'   }
%     {'R_RSC_ROI'   }
%     {'L_TPOJ2_ROI' }
%     {'R_TPOJ2_ROI' }
%     {'L_TPOJ3_ROI' }
%     {'R_TPOJ3_ROI' }
%     {'L_DVT_ROI'   }
%     {'R_DVT_ROI'   }
%     {'L_PGp_ROI'   }
%     {'R_PGp_ROI'   }
%     {'L_IP2_ROI'   }
%     {'R_IP2_ROI'   }
%     {'L_IP1_ROI'   }
%     {'R_IP1_ROI'   }
%     {'L_IP0_ROI'   }
%     {'R_IP0_ROI'   }
%     {'L_PFop_ROI'  }
%     {'R_PFop_ROI'  }
%     {'L_PF_ROI'    }
%     {'R_PF_ROI'    }
%     {'L_PFm_ROI'   }
%     {'R_PFm_ROI'   }
%     {'L_POS2_ROI'  }
%     {'R_POS2_ROI'  }
%     {'L_PGi_ROI'   }
%     {'R_PGi_ROI'   }
%     {'L_PGs_ROI'   }
%     {'R_PGs_ROI'   }
%     {'L_V6A_ROI'   }
%     {'R_V6A_ROI'   }
%     {'L_VMV1_ROI'  }
%     {'R_VMV1_ROI'  }
%     {'L_VMV3_ROI'  }
%     {'R_VMV3_ROI'  }
%     {'L_PHA2_ROI'  }
%     {'R_PHA2_ROI'  }
%     {'L_V4t_ROI'   }
%     {'R_V4t_ROI'   }
%     {'L_FST_ROI'   }
%     {'R_FST_ROI'   }
%     {'L_V3CD_ROI'  }
%     {'R_V3CD_ROI'  }
%     {'L_LO3_ROI'   }
%     {'R_LO3_ROI'   }
%     {'L_V7_ROI'    }
%     {'R_V7_ROI'    }
%     {'L_VMV2_ROI'  }
%     {'R_VMV2_ROI'  }
%     {'L_31pd_ROI'  }
%     {'R_31pd_ROI'  }
%     {'L_31a_ROI'   }
%     {'R_31a_ROI'   }
%     {'L_VVC_ROI'   }
%     {'R_VVC_ROI'   }
%     {'L_25_ROI'    }
%     {'R_25_ROI'    }
%     {'L_s32_ROI'   }
%     {'R_s32_ROI'   }
%     {'L_pOFC_ROI'  }
%     {'R_pOFC_ROI'  }
%     {'L_PoI1_ROI'  }
%     {'R_PoI1_ROI'  }
%     {'L_Ig_ROI'    }
%     {'R_Ig_ROI'    }
%     {'L_FOP5_ROI'  }
%     {'R_FOP5_ROI'  }
%     {'L_IPS1_ROI'  }
%     {'R_IPS1_ROI'  }
%     {'L_p10p_ROI'  }
%     {'R_p10p_ROI'  }
%     {'L_p47r_ROI'  }
%     {'R_p47r_ROI'  }
%     {'L_TGv_ROI'   }
%     {'R_TGv_ROI'   }
%     {'L_MBelt_ROI' }
%     {'R_MBelt_ROI' }
%     {'L_LBelt_ROI' }
%     {'R_LBelt_ROI' }
%     {'L_A4_ROI'    }
%     {'R_A4_ROI'    }
%     {'L_STSva_ROI' }
%     {'R_STSva_ROI' }
%     {'L_TE1m_ROI'  }
%     {'R_TE1m_ROI'  }
%     {'L_PI_ROI'    }
%     {'R_PI_ROI'    }
%     {'L_a32pr_ROI' }
%     {'R_a32pr_ROI' }
%     {'L_FFC_ROI'   }
%     {'R_FFC_ROI'   }
%     {'L_p24_ROI'   }
%     {'R_p24_ROI'   }
%     {'L_V3B_ROI'   }
%     {'R_V3B_ROI'   }
%     {'L_MST_ROI'   }
%     {'R_MST_ROI'   }
%     {'L_LO1_ROI'   }
%     {'R_LO1_ROI'   }
%     {'L_LO2_ROI'   }
%     {'R_LO2_ROI'   }
%     {'L_PIT_ROI'   }
%     {'R_PIT_ROI'   }
%     {'L_MT_ROI'    }
%     {'R_MT_ROI'    }
%     {'L_A1_ROI'    }
%     {'R_A1_ROI'    }
%     {'L_PSL_ROI'   }
%     {'R_PSL_ROI'   }
%     {'L_SFL_ROI'   }
%     {'R_SFL_ROI'   }
%     {'L_PCV_ROI'   }
%     {'R_PCV_ROI'   }
%     {'L_STV_ROI'   }
%     {'R_STV_ROI'   }
%     {'L_7Pm_ROI'   }
%     {'R_7Pm_ROI'   }
%     {'L_V6_ROI'    }
%     {'R_V6_ROI'    }
%     {'L_7m_ROI'    }
%     {'R_7m_ROI'    }
%     {'L_POS1_ROI'  }
%     {'R_POS1_ROI'  }
%     {'L_23d_ROI'   }
%     {'R_23d_ROI'   }
%     {'L_v23ab_ROI' }
%     {'R_v23ab_ROI' }
%     {'L_d23ab_ROI' }
%     {'R_d23ab_ROI' }
%     {'L_31pv_ROI'  }
%     {'R_31pv_ROI'  }
%     {'L_5m_ROI'    }
%     {'R_5m_ROI'    }
%     {'L_5mv_ROI'   }
%     {'R_5mv_ROI'   }
%     {'L_23c_ROI'   }
%     {'R_23c_ROI'   }
%     {'L_5L_ROI'    }
%     {'R_5L_ROI'    }
%     {'L_V2_ROI'    }
%     {'R_V2_ROI'    }
%     {'L_24dd_ROI'  }
%     {'R_24dd_ROI'  }
%     {'L_24dv_ROI'  }
%     {'R_24dv_ROI'  }
%     {'L_7AL_ROI'   }
%     {'R_7AL_ROI'   }
%     {'L_SCEF_ROI'  }
%     {'R_SCEF_ROI'  }
%     {'L_6ma_ROI'   }
%     {'R_6ma_ROI'   }
%     {'L_7Am_ROI'   }
%     {'R_7Am_ROI'   }
%     {'L_7PL_ROI'   }
%     {'R_7PL_ROI'   }
%     {'L_7PC_ROI'   }
%     {'R_7PC_ROI'   }
%     {'L_LIPv_ROI'  }
%     {'R_LIPv_ROI'  }
%     {'L_VIP_ROI'   }
%     {'R_VIP_ROI'   }
%     {'L_V3_ROI'    }
%     {'R_V3_ROI'    }
%     {'L_MIP_ROI'   }
%     {'R_MIP_ROI'   }
%     {'L_1_ROI'     }
%     {'R_1_ROI'     }
%     {'L_2_ROI'     }
%     {'R_2_ROI'     }
%     {'L_3a_ROI'    }
%     {'R_3a_ROI'    }
%     {'L_6d_ROI'    }
%     {'R_6d_ROI'    }
%     {'L_6mp_ROI'   }
%     {'R_6mp_ROI'   }
%     {'L_6v_ROI'    }
%     {'R_6v_ROI'    }
%     {'L_p24pr_ROI' }
%     {'R_p24pr_ROI' }
%     {'L_33pr_ROI'  }
%     {'R_33pr_ROI'  }
%     {'L_a24pr_ROI' }
%     {'R_a24pr_ROI' }
%     {'L_V4_ROI'    }
%     {'R_V4_ROI'    }
%     {'L_p32pr_ROI' }
%     {'R_p32pr_ROI' }
%     {'L_a24_ROI'   }
%     {'R_a24_ROI'   }
%     {'L_d32_ROI'   }
%     {'R_d32_ROI'   }
%     {'L_8BM_ROI'   }
%     {'R_8BM_ROI'   }
%     {'L_p32_ROI'   }
%     {'R_p32_ROI'   }
%     {'L_10r_ROI'   }
%     {'R_10r_ROI'   }
%     {'L_47m_ROI'   }
%     {'R_47m_ROI'   }
%     {'L_8Av_ROI'   }
%     {'R_8Av_ROI'   }
%     {'L_8Ad_ROI'   }
%     {'R_8Ad_ROI'   }
%     {'L_9m_ROI'    }
%     {'R_9m_ROI'    }
%     {'L_V8_ROI'    }
%     {'R_V8_ROI'    }
%     {'L_8BL_ROI'   }
%     {'R_8BL_ROI'   }
%     {'L_9p_ROI'    }
%     {'R_9p_ROI'    }
%     {'L_10d_ROI'   }
%     {'R_10d_ROI'   }
%     {'L_8C_ROI'    }
%     {'R_8C_ROI'    }
%     {'L_44_ROI'    }
%     {'R_44_ROI'    }
%     {'L_45_ROI'    }
%     {'R_45_ROI'    }
%     {'L_47l_ROI'   }
%     {'R_47l_ROI'   }
%     {'L_a47r_ROI'  }
%     {'R_a47r_ROI'  }
%     {'L_6r_ROI'    }
%     {'R_6r_ROI'    }
%     {'L_IFJa_ROI'  }
%     {'R_IFJa_ROI'  }
%     {'L_4_ROI'     }
%     {'R_4_ROI'     }
%     {'L_IFJp_ROI'  }
%     {'R_IFJp_ROI'  }
%     {'L_IFSp_ROI'  }
%     {'R_IFSp_ROI'  }
%     {'L_IFSa_ROI'  }
%     {'R_IFSa_ROI'  }
%     {'L_p9-46v_ROI'}
%     {'R_p9-46v_ROI'}
%     {'L_46_ROI'    }
%     {'R_46_ROI'    }
%     {'L_a9-46v_ROI'}
%     {'R_a9-46v_ROI'}
%     {'L_9-46d_ROI' }
%     {'R_9-46d_ROI' }
%     {'L_9a_ROI'    }
%     {'R_9a_ROI'    }
%     {'L_10v_ROI'   }
%     {'R_10v_ROI'   }
%     {'L_a10p_ROI'  }
%     {'R_a10p_ROI'  }
%     {'L_3b_ROI'    }
%     {'R_3b_ROI'    }
%     {'L_10pp_ROI'  }
%     {'R_10pp_ROI'  }
%     {'L_11l_ROI'   }
%     {'R_11l_ROI'   }
%     {'L_13l_ROI'   }
%     {'R_13l_ROI'   }
%     {'L_OFC_ROI'   }
%     {'R_OFC_ROI'   }
%     {'L_47s_ROI'   }
%     {'R_47s_ROI'   }
%     {'L_LIPd_ROI'  }
%     {'R_LIPd_ROI'  }
%     {'L_6a_ROI'    }
%     {'R_6a_ROI'    }
%     {'L_i6-8_ROI'  }
%     {'R_i6-8_ROI'  }
%     {'L_s6-8_ROI'  }
%     {'R_s6-8_ROI'  }
%     {'L_43_ROI'    }
%     {'R_43_ROI'    }

%
% region =
%     {'LO' }
%     {'RO' }
%     {'LC' }
%     {'RC' }
%     {'LC' }
%     {'RC' }
%     {'LC' }
%     {'RP' }
%     {'LC' }
%     {'RC' }
%     {'LT' }
%     {'RT' }
%     {'LT' }
%     {'RP' }
%     {'LP' }
%     {'RP' }
%     {'LT' }
%     {'RT' }
%     {'LT' }
%     {'RT' }
%     {'LF' }
%     {'RF' }
%     {'LT' }
%     {'RT' }
%     {'LC' }
%     {'RF' }
%     {'LT' }
%     {'RT' }
%     {'LPF'}
%     {'RPF'}
%     {'LT' }
%     {'RT' }
%     {'LC' }
%     {'RC' }
%     {'LT' }
%     {'RT' }
%     {'LC' }
%     {'RC' }
%     {'LP' }
%     {'RP' }
%     {'LP' }
%     {'RP' }
%     {'UU' }
%     {'UU' }
%     {'UU' }
%     {'UU' }
%     {'LC' }
%     {'RC' }
%     {'UU' }
%     {'UU' }
%     {'LO' }
%     {'RO' }
%     {'LT' }
%     {'RT' }
%     {'LT' }
%     {'RT' }
%     {'LT' }
%     {'RT' }
%     {'LT' }
%     {'RT' }
%     {'LT' }
%     {'RT' }
%     {'LT' }
%     {'RT' }
%     {'LT' }
%     {'RT' }
%     {'LT' }
%     {'RT' }
%     {'LP' }
%     {'RP' }
%     {'LT' }
%     {'RT' }
%     {'LT' }
%     {'RT' }
%     {'LT' }
%     {'RT' }
%     {'LT' }
%     {'RT' }
%     {'LT' }
%     {'RT' }
%     {'LT' }
%     {'RT' }
%     {'LT' }
%     {'RT' }
%     {'LT' }
%     {'RT' }
%     {'LT' }
%     {'RT' }
%     {'LT' }
%     {'RT' }
%     {'LL' }
%     {'RL' }
%     {'LP' }
%     {'RT' }
%     {'LP' }
%     {'RP' }
%     {'LP' }
%     {'RP' }
%     {'LP' }
%     {'RP' }
%     {'LP' }
%     {'RP' }
%     {'LP' }
%     {'RP' }
%     {'LP' }
%     {'RP' }
%     {'LP' }
%     {'RP' }
%     {'LP' }
%     {'RP' }
%     {'LP' }
%     {'RP' }
%     {'LP' }
%     {'RP' }
%     {'LP' }
%     {'RP' }
%     {'LP' }
%     {'RP' }
%     {'LP' }
%     {'RP' }
%     {'LO' }
%     {'RO' }
%     {'LT' }
%     {'RT' }
%     {'LT' }
%     {'RT' }
%     {'LO' }
%     {'RO' }
%     {'LO' }
%     {'RT' }
%     {'LO' }
%     {'RO' }
%     {'LO' }
%     {'RP' }
%     {'LP' }
%     {'RP' }
%     {'LO' }
%     {'RO' }
%     {'LP' }
%     {'RP' }
%     {'LP' }
%     {'RP' }
%     {'LT' }
%     {'RT' }
%     {'LP' }
%     {'RP' }
%     {'LPF'}
%     {'RPF'}
%     {'LPF'}
%     {'RPF'}
%     {'LT' }
%     {'RT' }
%     {'LT' }
%     {'RT' }
%     {'LF' }
%     {'RF' }
%     {'LP' }
%     {'RP' }
%     {'LF' }
%     {'RF' }
%     {'LF' }
%     {'RF' }
%     {'LT' }
%     {'RT' }
%     {'LT' }
%     {'RT' }
%     {'LT' }
%     {'RT' }
%     {'LT' }
%     {'RT' }
%     {'LT' }
%     {'RT' }
%     {'LT' }
%     {'RT' }
%     {'LT' }
%     {'RT' }
%     {'LF' }
%     {'RF' }
%     {'LT' }
%     {'RT' }
%     {'LL' }
%     {'RL' }
%     {'LP' }
%     {'RP' }
%     {'LO' }
%     {'RO' }
%     {'LO' }
%     {'RO' }
%     {'LO' }
%     {'RO' }
%     {'LO' }
%     {'RO' }
%     {'LO' }
%     {'RO' }
%     {'LT' }
%     {'RT' }
%     {'LP' }
%     {'RP' }
%     {'LF' }
%     {'RF' }
%     {'LP' }
%     {'RP' }
%     {'LT' }
%     {'RP' }
%     {'LP' }
%     {'RP' }
%     {'LP' }
%     {'RP' }
%     {'LP' }
%     {'RP' }
%     {'LP' }
%     {'RP' }
%     {'LL' }
%     {'RL' }
%     {'LP' }
%     {'RL' }
%     {'LL' }
%     {'RL' }
%     {'LP' }
%     {'RP' }
%     {'LC' }
%     {'RC' }
%     {'LC' }
%     {'RC' }
%     {'LL' }
%     {'RP' }
%     {'LP' }
%     {'RP' }
%     {'LO' }
%     {'RO' }
%     {'LF' }
%     {'RC' }
%     {'LF' }
%     {'RF' }
%     {'LP' }
%     {'RP' }
%     {'LF' }
%     {'RF' }
%     {'LF' }
%     {'RF' }
%     {'LP' }
%     {'RP' }
%     {'LP' }
%     {'RP' }
%     {'LP' }
%     {'RP' }
%     {'LP' }
%     {'RP' }
%     {'LP' }
%     {'RP' }
%     {'LO' }
%     {'RO' }
%     {'LP' }
%     {'RP' }
%     {'LC' }
%     {'RC' }
%     {'LC' }
%     {'RC' }
%     {'LC' }
%     {'RC' }
%     {'LC' }
%     {'RC' }
%     {'LC' }
%     {'RC' }
%     {'LC' }
%     {'RC' }
%     {'LL' }
%     {'RL' }
%     {'LL' }
%     {'RL' }
%     {'LL' }
%     {'RL' }
%     {'LO' }
%     {'RO' }
%     {'LF' }
%     {'RF' }
%     {'LL' }
%     {'RL' }
%     {'LF' }
%     {'RF' }
%     {'LF' }
%     {'RF' }
%     {'LPF'}
%     {'RPF'}
%     {'LPF'}
%     {'RPF'}
%     {'LPF'}
%     {'RPF'}
%     {'LF' }
%     {'RF' }
%     {'LF' }
%     {'RF' }
%     {'LF' }
%     {'RF' }
%     {'LT' }
%     {'RT' }
%     {'LF' }
%     {'RF' }
%     {'LF' }
%     {'RF' }
%     {'LF' }
%     {'RF' }
%     {'LF' }
%     {'RF' }
%     {'LF' }
%     {'RF' }
%     {'LF' }
%     {'RF' }
%     {'LPF'}
%     {'RPF'}
%     {'LPF'}
%     {'RF' }
%     {'LC' }
%     {'RF' }
%     {'LF' }
%     {'RF' }
%     {'LC' }
%     {'RC' }
%     {'LF' }
%     {'RF' }
%     {'LF' }
%     {'RF' }
%     {'LF' }
%     {'RF' }
%     {'LF' }
%     {'RF' }
%     {'LF' }
%     {'RF' }
%     {'LF' }
%     {'RF' }
%     {'LF' }
%     {'RF' }
%     {'LF' }
%     {'RF' }
%     {'LP' }
%     {'RP' }
%     {'LF' }
%     {'RF' }
%     {'LC' }
%     {'RC' }
%     {'LPF'}
%     {'RPF'}
%     {'LPF'}
%     {'RPF'}
%     {'LPF'}
%     {'RPF'}
%     {'LPF'}
%     {'RPF'}
%     {'LPF'}
%     {'RPF'}
%     {'LP' }
%     {'RP' }
%     {'LF' }
%     {'RF' }
%     {'LF' }
%     {'RF' }
%     {'LF' }
%     {'RF' }
%     {'LC' }
%     {'RC' }
