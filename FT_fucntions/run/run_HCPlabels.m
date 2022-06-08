% groups = {{'Primary Visual Cortex (V1)','V1'}, ...
%     {'Early Visual Cortex','V2', 'V3', 'V4'}, ...
%     {'Dorsal Stream Visual Cortex','V3A', 'V3B', 'V6', 'V6A', 'V7', 'IPS1'}, ...
%     {'Ventral Stream Visual Cortex','V8', 'VVC', 'PIT', 'FFC', 'VMV1', 'VMV2', 'VMV3'}, ...
%     {'MT+ Complex and Neighboring Visual Areas', 'V3CD', 'LO1', 'LO2', 'LO3', 'V4t', 'FST', 'MT', 'MST', 'PH'},...
%     {'Somatosensory and Motor Cortex','4', '3a', '3b', '1', '2'},...
%     {'Paracentral Lobular and Mid Cingulate Cortex','24dd', '24dv', '6mp', '6ma', 'SCEF', '5m', '5L', '5mv'},...
%     {'Premotor Cortex','55b', '6d', '6a', 'FEF', '6v', '6r', 'PEF'},...
%     {'Posterior Opercular Cortex','43', 'FOP1', 'OP4', 'OP1', 'OP2-3', 'PFcm'},...
%     {'Early Auditory Cortex','A1', 'LBelt', 'MBelt', 'PBelt', 'RI'},...
%     {'Auditory Association Cortex','A4', 'A5', 'STSdp', 'STSda', 'STSvp', 'STSva', 'STGa', 'TA2'},...
%     {'Insular and Frontal Opercular Cortex','52', 'PI', 'Ig', 'PoI1', 'PoI2', 'FOP2', 'FOP3', ...
%     'MI', 'AVI', 'AAIC', 'Pir', 'FOP4', 'FOP5'},...
%     {'Medial Temporal Cortex','H', 'PreS', 'EC', 'PeEc', 'PHA1', 'PHA2', 'PHA3'},...
%     {'Lateral Temporal Cortex','PHT', 'TE1p', 'TE1m', 'TE1a', 'TE2p', 'TE2a', ...
%     'TGv', 'TGd', 'TF'},...
%     {'Temporo-Parieto-Occipital Junction','TPOJ1', 'TPOJ2', 'TPOJ3', 'STV', 'PSL'},...
%     {'Superior Parietal Cortex','LIPv', 'LIPd', 'VIP', 'AIP', 'MIP', ...
%     '7PC', '7AL', '7Am', '7PL', '7Pm'},...
%     {'Inferior Parietal Cortex','PGp', 'PGs', 'PGi', 'PFm', 'PF', 'PFt', 'PFop', ...
%     'IP0', 'IP1', 'IP2'},...
%     {'Posterior Cingulate Cortex','DVT', 'ProS', 'POS1', 'POS2', 'RSC', 'v23ab', 'd23ab', ...
%     '31pv', '31pd', '31a', '23d', '23c', 'PCV', '7m'},...
%     {'Anterior Cingulate and Medial Prefrontal Cortex','33pr', 'p24pr', 'a24pr', 'p24', 'a24', 'p32pr', 'a32pr', 'd32', ...
%     'p32', 's32', '8BM', '9m', '10v', '10r', '25',},...
%     {'Orbital and Polar Frontal Cortex','47s', '47m', 'a47r', '11l', '13l',...
%     'a10p', 'p10p', '10pp', '10d', 'OFC', 'pOFC',},...
%     {'Inferior Frontal Cortex','44', '45', 'IFJp', 'IFJa', 'IFSp', 'IFSa', '47l', 'p47r',},...
%     {'DorsoLateral Prefrontal Cortex','8C', '8Av', 'i6-8', 's6-8', 'SFL', '8BL', '9p', '9a', '8Ad',...
%     'p9-46v', 'a9-46v', '46', '9-46d',}, ...
%     {'???','???'}};

% groups = {'?-lh';'?-rh';'10d_ROI-lh';'10d_ROI-rh';'10pp_ROI-lh';'10pp_ROI-rh';'10r_ROI-lh';'10r_ROI-rh';'10v_ROI-lh';'10v_ROI-rh';'11l_ROI-lh';'11l_ROI-rh';'13l_ROI-lh';'13l_ROI-rh';'1_ROI-lh';'1_ROI-rh';'23c_ROI-lh';'23c_ROI-rh';'23d_ROI-lh';'23d_ROI-rh';'24dd_ROI-lh';'24dd_ROI-rh';'24dv_ROI-lh';'24dv_ROI-rh';'25_ROI-lh';'25_ROI-rh';'2_ROI-lh';'2_ROI-rh';'31a_ROI-lh';'31a_ROI-rh';'31pd_ROI-lh';'31pd_ROI-rh';'31pv_ROI-lh';'31pv_ROI-rh';'33pr_ROI-lh';'33pr_ROI-rh';'3a_ROI-lh';'3a_ROI-rh';'3b_ROI-lh';'3b_ROI-rh';'43_ROI-lh';'43_ROI-rh';'44_ROI-lh';'44_ROI-rh';'45_ROI-lh';'45_ROI-rh';'46_ROI-lh';'46_ROI-rh';'47l_ROI-lh';'47l_ROI-rh';'47m_ROI-lh';'47m_ROI-rh';'47s_ROI-lh';'47s_ROI-rh';'4_ROI-lh';'4_ROI-rh';'52_ROI-lh';'52_ROI-rh';'55b_ROI-lh';'55b_ROI-rh';'5L_ROI-lh';'5L_ROI-rh';'5m_ROI-lh';'5m_ROI-rh';'5mv_ROI-lh';'5mv_ROI-rh';'6a_ROI-lh';'6a_ROI-rh';'6d_ROI-lh';'6d_ROI-rh';'6ma_ROI-lh';'6ma_ROI-rh';'6mp_ROI-lh';'6mp_ROI-rh';'6r_ROI-lh';'6r_ROI-rh';'6v_ROI-lh';'6v_ROI-rh';'7AL_ROI-lh';'7AL_ROI-rh';'7Am_ROI-lh';'7Am_ROI-rh';'7PC_ROI-lh';'7PC_ROI-rh';'7PL_ROI-lh';'7PL_ROI-rh';'7Pm_ROI-lh';'7Pm_ROI-rh';'7m_ROI-lh';'7m_ROI-rh';'8Ad_ROI-lh';'8Ad_ROI-rh';'8Av_ROI-lh';'8Av_ROI-rh';'8BL_ROI-lh';'8BL_ROI-rh';'8BM_ROI-lh';'8BM_ROI-rh';'8C_ROI-lh';'8C_ROI-rh';'9-46d_ROI-lh';'9-46d_ROI-rh';'9a_ROI-lh';'9a_ROI-rh';'9m_ROI-lh';'9m_ROI-rh';'9p_ROI-lh';'9p_ROI-rh';'A1_ROI-lh';'A1_ROI-rh';'A4_ROI-lh';'A4_ROI-rh';'A5_ROI-lh';'A5_ROI-rh';'AAIC_ROI-lh';'AAIC_ROI-rh';'AIP_ROI-lh';'AIP_ROI-rh';'AVI_ROI-lh';'AVI_ROI-rh';'DVT_ROI-lh';'DVT_ROI-rh';'EC_ROI-lh';'EC_ROI-rh';'FEF_ROI-lh';'FEF_ROI-rh';'FFC_ROI-lh';'FFC_ROI-rh';'FOP1_ROI-lh';'FOP1_ROI-rh';'FOP2_ROI-lh';'FOP2_ROI-rh';'FOP3_ROI-lh';'FOP3_ROI-rh';'FOP4_ROI-lh';'FOP4_ROI-rh';'FOP5_ROI-lh';'FOP5_ROI-rh';'FST_ROI-lh';'FST_ROI-rh';'H_ROI-lh';'H_ROI-rh';'IFJa_ROI-lh';'IFJa_ROI-rh';'IFJp_ROI-lh';'IFJp_ROI-rh';'IFSa_ROI-lh';'IFSa_ROI-rh';'IFSp_ROI-lh';'IFSp_ROI-rh';'IP0_ROI-lh';'IP0_ROI-rh';'IP1_ROI-lh';'IP1_ROI-rh';'IP2_ROI-lh';'IP2_ROI-rh';'IPS1_ROI-lh';'IPS1_ROI-rh';'Ig_ROI-lh';'Ig_ROI-rh';'LBelt_ROI-lh';'LBelt_ROI-rh';'LIPd_ROI-lh';'LIPd_ROI-rh';'LIPv_ROI-lh';'LIPv_ROI-rh';'LO1_ROI-lh';'LO1_ROI-rh';'LO2_ROI-lh';'LO2_ROI-rh';'LO3_ROI-lh';'LO3_ROI-rh';'MBelt_ROI-lh';'MBelt_ROI-rh';'MIP_ROI-lh';'MIP_ROI-rh';'MI_ROI-lh';'MI_ROI-rh';'MST_ROI-lh';'MST_ROI-rh';'MT_ROI-lh';'MT_ROI-rh';'OFC_ROI-lh';'OFC_ROI-rh';'OP1_ROI-lh';'OP1_ROI-rh';'OP2-3_ROI-lh';'OP2-3_ROI-rh';'OP4_ROI-lh';'OP4_ROI-rh';'PBelt_ROI-lh';'PBelt_ROI-rh';'PCV_ROI-lh';'PCV_ROI-rh';'PEF_ROI-lh';'PEF_ROI-rh';'PF_ROI-lh';'PF_ROI-rh';'PFcm_ROI-lh';'PFcm_ROI-rh';'PFm_ROI-lh';'PFm_ROI-rh';'PFop_ROI-lh';'PFop_ROI-rh';'PFt_ROI-lh';'PFt_ROI-rh';'PGi_ROI-lh';'PGi_ROI-rh';'PGp_ROI-lh';'PGp_ROI-rh';'PGs_ROI-lh';'PGs_ROI-rh';'PHA1_ROI-lh';'PHA1_ROI-rh';'PHA2_ROI-lh';'PHA2_ROI-rh';'PHA3_ROI-lh';'PHA3_ROI-rh';'PHT_ROI-lh';'PHT_ROI-rh';'PH_ROI-lh';'PH_ROI-rh';'PIT_ROI-lh';'PIT_ROI-rh';'PI_ROI-lh';'PI_ROI-rh';'POS1_ROI-lh';'POS1_ROI-rh';'POS2_ROI-lh';'POS2_ROI-rh';'PSL_ROI-lh';'PSL_ROI-rh';'PeEc_ROI-lh';'PeEc_ROI-rh';'Pir_ROI-lh';'Pir_ROI-rh';'PoI1_ROI-lh';'PoI1_ROI-rh';'PoI2_ROI-lh';'PoI2_ROI-rh';'PreS_ROI-lh';'PreS_ROI-rh';'ProS_ROI-lh';'ProS_ROI-rh';'RI_ROI-lh';'RI_ROI-rh';'RSC_ROI-lh';'RSC_ROI-rh';'SCEF_ROI-lh';'SCEF_ROI-rh';'SFL_ROI-lh';'SFL_ROI-rh';'STGa_ROI-lh';'STGa_ROI-rh';'STSda_ROI-lh';'STSda_ROI-rh';'STSdp_ROI-lh';'STSdp_ROI-rh';'STSva_ROI-lh';'STSva_ROI-rh';'STSvp_ROI-lh';'STSvp_ROI-rh';'STV_ROI-lh';'STV_ROI-rh';'TA2_ROI-lh';'TA2_ROI-rh';'TE1a_ROI-lh';'TE1a_ROI-rh';'TE1m_ROI-lh';'TE1m_ROI-rh';'TE1p_ROI-lh';'TE1p_ROI-rh';'TE2a_ROI-lh';'TE2a_ROI-rh';'TE2p_ROI-lh';'TE2p_ROI-rh';'TF_ROI-lh';'TF_ROI-rh';'TGd_ROI-lh';'TGd_ROI-rh';'TGv_ROI-lh';'TGv_ROI-rh';'TPOJ1_ROI-lh';'TPOJ1_ROI-rh';'TPOJ2_ROI-lh';'TPOJ2_ROI-rh';'TPOJ3_ROI-lh';'TPOJ3_ROI-rh';'V1_ROI-lh';'V1_ROI-rh';'V2_ROI-lh';'V2_ROI-rh';'V3A_ROI-lh';'V3A_ROI-rh';'V3B_ROI-lh';'V3B_ROI-rh';'V3CD_ROI-lh';'V3CD_ROI-rh';'V3_ROI-lh';'V3_ROI-rh';'V4_ROI-lh';'V4_ROI-rh';'V4t_ROI-lh';'V4t_ROI-rh';'V6A_ROI-lh';'V6A_ROI-rh';'V6_ROI-lh';'V6_ROI-rh';'V7_ROI-lh';'V7_ROI-rh';'V8_ROI-lh';'V8_ROI-rh';'VIP_ROI-lh';'VIP_ROI-rh';'VMV1_ROI-lh';'VMV1_ROI-rh';'VMV2_ROI-lh';'VMV2_ROI-rh';'VMV3_ROI-lh';'VMV3_ROI-rh';'VVC_ROI-lh';'VVC_ROI-rh';'a10p_ROI-lh';'a10p_ROI-rh';'a24_ROI-lh';'a24_ROI-rh';'a24pr_ROI-lh';'a24pr_ROI-rh';'a32pr_ROI-lh';'a32pr_ROI-rh';'a47r_ROI-lh';'a47r_ROI-rh';'a9-46v_ROI-lh';'a9-46v_ROI-rh';'d23ab_ROI-lh';'d23ab_ROI-rh';'d32_ROI-lh';'d32_ROI-rh';'i6-8_ROI-lh';'i6-8_ROI-rh';'p10p_ROI-lh';'p10p_ROI-rh';'p24_ROI-lh';'p24_ROI-rh';'p24pr_ROI-lh';'p24pr_ROI-rh';'p32_ROI-lh';'p32_ROI-rh';'p32pr_ROI-lh';'p32pr_ROI-rh';'p47r_ROI-lh';'p47r_ROI-rh';'p9-46v_ROI-lh';'p9-46v_ROI-rh';'pOFC_ROI-lh';'pOFC_ROI-rh';'s32_ROI-lh';'s32_ROI-rh';'s6-8_ROI-lh';'s6-8_ROI-rh';'v23ab_ROI-lh';'v23ab_ROI-rh'};

groups = {'???_ROI';         
    'V1_ROI';    
    'MST_ROI';   
    'V6_ROI';    
    'V2_ROI';    
    'V3_ROI';    
    'V4_ROI';    
    'V8_ROI';    
    '4_ROI';     
    '3b_ROI';    
    'FEF_ROI';   
    'PEF_ROI';   
    '55b_ROI';   
    'V3A_ROI';   
    'RSC_ROI';   
    'POS2_ROI';  
    'V7_ROI';    
    'IPS1_ROI';  
    'FFC_ROI';   
    'V3B_ROI';   
    'LO1_ROI';   
    'LO2_ROI';   
    'PIT_ROI';   
    'MT_ROI';    
    'A1_ROI';    
    'PSL_ROI';   
    'SFL_ROI';   
    'PCV_ROI';   
    'STV_ROI';   
    '7Pm_ROI';   
    '7m_ROI';    
    'POS1_ROI';  
    '23d_ROI';   
    'v23ab_ROI'; 
    'd23ab_ROI'; 
    '31pv_ROI';  
    '5m_ROI';    
    '5mv_ROI';   
    '23c_ROI';   
    '5L_ROI';    
    '24dd_ROI';  
    '24dv_ROI';  
    '7AL_ROI';   
    'SCEF_ROI';  
    '6ma_ROI';   
    '7Am_ROI';   
    '7PL_ROI';   
    '7PC_ROI';   
    'LIPv_ROI';  
    'VIP_ROI';   
    'MIP_ROI';   
    '1_ROI';     
    '2_ROI';     
    '3a_ROI';    
    '6d_ROI';    
    '6mp_ROI';   
    '6v_ROI';    
    'p24pr_ROI'; 
    '33pr_ROI';  
    'a24pr_ROI'; 
    'p32pr_ROI'; 
    'a24_ROI';   
    'd32_ROI';   
    '8BM_ROI';   
    'p32_ROI';   
    '10r_ROI';   
    '47m_ROI';   
    '8Av_ROI';   
    '8Ad_ROI';   
    '9m_ROI';    
    '8BL_ROI';   
    '9p_ROI';    
    '10d_ROI';   
    '8C_ROI';    
    '44_ROI';    
    '45_ROI';    
    '47l_ROI';   
    'a47r_ROI';  
    '6r_ROI';    
    'IFJa_ROI';  
    'IFJp_ROI';  
    'IFSp_ROI';  
    'IFSa_ROI';  
    'p9-46v_ROI';
    '46_ROI';    
    'a9-46v_ROI';
    '9-46d_ROI'; 
    '9a_ROI';    
    '10v_ROI';   
    'a10p_ROI';  
    '10pp_ROI';  
    '11l_ROI';   
    '13l_ROI';   
    'OFC_ROI';   
    '47s_ROI';   
    'LIPd_ROI';  
    '6a_ROI';    
    'i6-8_ROI';  
    's6-8_ROI';  
    '43_ROI';    
    'OP4_ROI';   
    'OP1_ROI';   
    'OP2-3_ROI'; 
    '52_ROI';    
    'RI_ROI';    
    'PFcm_ROI';  
    'PoI2_ROI';  
    'TA2_ROI';   
    'FOP4_ROI';  
    'MI_ROI';    
    'Pir_ROI';   
    'AVI_ROI';   
    'AAIC_ROI';  
    'FOP1_ROI';  
    'FOP3_ROI';  
    'FOP2_ROI';  
    'PFt_ROI';   
    'AIP_ROI';   
    'EC_ROI';    
    'PreS_ROI';  
    'H_ROI';     
    'ProS_ROI';  
    'PeEc_ROI';  
    'STGa_ROI';  
    'PBelt_ROI'; 
    'A5_ROI';    
    'PHA1_ROI';  
    'PHA3_ROI';  
    'STSda_ROI'; 
    'STSdp_ROI'; 
    'STSvp_ROI'; 
    'TGd_ROI';   
    'TE1a_ROI';  
    'TE1p_ROI';  
    'TE2a_ROI';  
    'TF_ROI';    
    'TE2p_ROI';  
    'PHT_ROI';   
    'PH_ROI';    
    'TPOJ1_ROI'; 
    'TPOJ2_ROI'; 
    'TPOJ3_ROI'; 
    'DVT_ROI';   
    'PGp_ROI';   
    'IP2_ROI';   
    'IP1_ROI';   
    'IP0_ROI';   
    'PFop_ROI';  
    'PF_ROI';    
    'PFm_ROI';   
    'PGi_ROI';   
    'PGs_ROI';   
    'V6A_ROI';   
    'VMV1_ROI';  
    'VMV3_ROI';  
    'PHA2_ROI';  
    'V4t_ROI';   
    'FST_ROI';   
    'V3CD_ROI';  
    'LO3_ROI';   
    'VMV2_ROI';  
    '31pd_ROI';  
    '31a_ROI';   
    'VVC_ROI';   
    '25_ROI';    
    's32_ROI';   
    'pOFC_ROI';  
    'PoI1_ROI';  
    'Ig_ROI';    
    'FOP5_ROI';  
    'p10p_ROI';  
    'p47r_ROI';  
    'TGv_ROI';   
    'MBelt_ROI'; 
    'LBelt_ROI'; 
    'A4_ROI';    
    'STSva_ROI'; 
    'TE1m_ROI';  
    'PI_ROI';    
    'a32pr_ROI'; 
    'p24_ROI'};   
% {'???_ROI';'L_V1_ROI_ROI';'L_MST_ROI_ROI';'L_V6_ROI_ROI';'L_V2_ROI_ROI';'L_V3_ROI_ROI';'L_V4_ROI_ROI';'L_V8_ROI_ROI';'L_4_ROI_ROI';'L_3b_ROI_ROI';'L_FEF_ROI_ROI';'L_PEF_ROI_ROI';'L_55b_ROI_ROI';'L_V3A_ROI_ROI';'L_RSC_ROI_ROI';'L_POS2_ROI_ROI';'L_V7_ROI_ROI';'L_IPS1_ROI_ROI';'L_FFC_ROI_ROI';'L_V3B_ROI_ROI';'L_LO1_ROI_ROI';'L_LO2_ROI_ROI';'L_PIT_ROI_ROI';'L_MT_ROI_ROI';'L_A1_ROI_ROI';'L_PSL_ROI_ROI';'L_SFL_ROI_ROI';'L_PCV_ROI_ROI';'L_STV_ROI_ROI';'L_7Pm_ROI_ROI';'L_7m_ROI_ROI';'L_POS1_ROI_ROI';'L_23d_ROI_ROI';'L_v23ab_ROI_ROI';'L_d23ab_ROI_ROI';'L_31pv_ROI_ROI';'L_5m_ROI_ROI';'L_5mv_ROI_ROI';'L_23c_ROI_ROI';'L_5L_ROI_ROI';'L_24dd_ROI_ROI';'L_24dv_ROI_ROI';'L_7AL_ROI_ROI';'L_SCEF_ROI_ROI';'L_6ma_ROI_ROI';'L_7Am_ROI_ROI';'L_7PL_ROI_ROI';'L_7PC_ROI_ROI';'L_LIPv_ROI_ROI';'L_VIP_ROI_ROI';'L_MIP_ROI_ROI';'L_1_ROI_ROI';'L_2_ROI_ROI';'L_3a_ROI_ROI';'L_6d_ROI_ROI';'L_6mp_ROI_ROI';'L_6v_ROI_ROI';'L_p24pr_ROI_ROI';'L_33pr_ROI_ROI';'L_a24pr_ROI_ROI';'L_p32pr_ROI_ROI';'L_a24_ROI_ROI';'L_d32_ROI_ROI';'L_8BM_ROI_ROI';'L_p32_ROI_ROI';'L_10r_ROI_ROI';'L_47m_ROI_ROI';'L_8Av_ROI_ROI';'L_8Ad_ROI_ROI';'L_9m_ROI_ROI';'L_8BL_ROI_ROI';'L_9p_ROI_ROI';'L_10d_ROI_ROI';'L_8C_ROI_ROI';'L_44_ROI_ROI';'L_45_ROI_ROI';'L_47l_ROI_ROI';'L_a47r_ROI_ROI';'L_6r_ROI_ROI';'L_IFJa_ROI_ROI';'L_IFJp_ROI_ROI';'L_IFSp_ROI_ROI';'L_IFSa_ROI_ROI';'L_p9-46v_ROI_ROI';'L_46_ROI_ROI';'L_a9-46v_ROI_ROI';'L_9-46d_ROI_ROI';'L_9a_ROI_ROI';'L_10v_ROI_ROI';'L_a10p_ROI_ROI';'L_10pp_ROI_ROI';'L_11l_ROI_ROI';'L_13l_ROI_ROI';'L_OFC_ROI_ROI';'L_47s_ROI_ROI';'L_LIPd_ROI_ROI';'L_6a_ROI_ROI';'L_i6-8_ROI_ROI';'L_s6-8_ROI_ROI';'L_43_ROI_ROI';'L_OP4_ROI_ROI';'L_OP1_ROI_ROI';'L_OP2-3_ROI_ROI';'L_52_ROI_ROI';'L_RI_ROI_ROI';'L_PFcm_ROI_ROI';'L_PoI2_ROI_ROI';'L_TA2_ROI_ROI';'L_FOP4_ROI_ROI';'L_MI_ROI_ROI';'L_Pir_ROI_ROI';'L_AVI_ROI_ROI';'L_AAIC_ROI_ROI';'L_FOP1_ROI_ROI';'L_FOP3_ROI_ROI';'L_FOP2_ROI_ROI';'L_PFt_ROI_ROI';'L_AIP_ROI_ROI';'L_EC_ROI_ROI';'L_PreS_ROI_ROI';'L_H_ROI_ROI';'L_ProS_ROI_ROI';'L_PeEc_ROI_ROI';'L_STGa_ROI_ROI';'L_PBelt_ROI_ROI';'L_A5_ROI_ROI';'L_PHA1_ROI_ROI';'L_PHA3_ROI_ROI';'L_STSda_ROI_ROI';'L_STSdp_ROI_ROI';'L_STSvp_ROI_ROI';'L_TGd_ROI_ROI';'L_TE1a_ROI_ROI';'L_TE1p_ROI_ROI';'L_TE2a_ROI_ROI';'L_TF_ROI_ROI';'L_TE2p_ROI_ROI';'L_PHT_ROI_ROI';'L_PH_ROI_ROI';'L_TPOJ1_ROI_ROI';'L_TPOJ2_ROI_ROI';'L_TPOJ3_ROI_ROI';'L_DVT_ROI_ROI';'L_PGp_ROI_ROI';'L_IP2_ROI_ROI';'L_IP1_ROI_ROI';'L_IP0_ROI_ROI';'L_PFop_ROI_ROI';'L_PF_ROI_ROI';'L_PFm_ROI_ROI';'L_PGi_ROI_ROI';'L_PGs_ROI_ROI';'L_V6A_ROI_ROI';'L_VMV1_ROI_ROI';'L_VMV3_ROI_ROI';'L_PHA2_ROI_ROI';'L_V4t_ROI_ROI';'L_FST_ROI_ROI';'L_V3CD_ROI_ROI';'L_LO3_ROI_ROI';'L_VMV2_ROI_ROI';'L_31pd_ROI_ROI';'L_31a_ROI_ROI';'L_VVC_ROI_ROI';'L_25_ROI_ROI';'L_s32_ROI_ROI';'L_pOFC_ROI_ROI';'L_PoI1_ROI_ROI';'L_Ig_ROI_ROI';'L_FOP5_ROI_ROI';'L_p10p_ROI_ROI';'L_p47r_ROI_ROI';'L_TGv_ROI_ROI';'L_MBelt_ROI_ROI';'L_LBelt_ROI_ROI';'L_A4_ROI_ROI';'L_STSva_ROI_ROI';'L_TE1m_ROI_ROI';'L_PI_ROI_ROI';'L_a32pr_ROI_ROI';'L_p24_ROI'}
%%
groups_labels = groups;

% groups_labels = {{'Prim Visual C. (V1)'}, ...
%     {'Early Vis C.'}, ...
%     {'Dors Stream Vis C.'}, ...
%     {'Vent Stream Vis C.'}, ...
%     {'MT+ Complex and Neighboring Vis Areas'},...
%     {'Somatosensory and Motor C.'},...
%     {'Paracentral Lobular and Mid Cin C.',},...
%     {'Premotor C.'},...
%     {'Post Opercular C.'},...
%     {'Early Aud C.'},...
%     {'Aud Association C.'},...
%     {'Insular and Front Oper C.'},...
%     {'Med Temp C.'},...
%     {'Lateral Temp C.'},...
%     {'Temporo-Parieto-Occip Junction'},...
%     {'Inf Pari C.'},...
%     {'Post Cin C.'},...
%     {'Ant Cin and Med PreFront C.'},...
%     {'Orbital and Polar Front C.'},...
%     {'Inf Front C.'},...
%     {'DorsoLateral PreFront C.'}, ...
%     {'???','???'}};

