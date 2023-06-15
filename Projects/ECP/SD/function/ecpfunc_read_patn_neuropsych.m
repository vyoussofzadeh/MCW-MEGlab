function patn_neuropsych_data = ecpfunc_read_patn_neuropsych()

%% Setup the Import Options
opts = delimitedTextImportOptions("NumVariables", 75);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["SUBNO", "Group", "Cogntv_Group", "Sex", "Age", "Ethnic", "Race", "Degree", "EducationYRS", "DomntHand", "EHQ", "DegreeMother", "EDyrsMother", "DegreeFather", "EDyrsFather", "FamHistEpi", "AgeFirstSz", "AgeRecurrSzs", "AgeAEDstart", "AEDcount", "TLESide", "SeizureTypesList", "Aura_OnsetAge", "Aura_MonthFreq", "SP_OnsetAge", "SP_MonthFreq", "CP_OnsetAge", "CP_MonthFreq", "GTC_OnsetAge", "GTC_MonthFreq", "RecentCP_days", "RecentGTC_days", "EEG_InterIctalActv", "EEG_InterictalLOC", "EEG_IctalActv", "EEG_IctalLOC", "EEGLateralization", "WASI_BlckR_ZScore", "WASI_VocR_ZScore", "JOLO_R_ZScore", "GrovPegD_R_ZScore", "GrovPegND_R_ZScore", "COWA_R_ZScore", "SemntFl_R_ZScore", "RAVLT_TotalR_ZScore", "RAVLT_DelayRecalR_ZScore", "BNT_R_ZScore", "DCCS_raw_ZScore", "DCCS_computed_ZScore", "FLANKv21_raw_ZScore", "FLANKv21_computed_ZScore", "WMEM_raw_ZScore", "ODOR_raw_ZScore", "ORALR_theta_ZScore", "PSPEED_raw_ZScore", "PSPEED_computed_ZScore", "SQMEM_raw_ZScore", "SQMEM_computed_ZScore", "SQMEM_theta_ZScore", "VOCAB_theta_ZScore", "GPD_R_Z_inv", "GPND_R_Z_inv", "BFRT_Zscore", "Language", "EF_Speed", "Memory", "Visuospatial", "Motor_Speed", "Memory2", "Languagenew", "Visuospatialnew", "EF_VMSpeednew", "EF_Verbalnew", "Memorynew", "Motornew"];
opts.VariableTypes = ["double", "categorical", "double", "categorical", "double", "categorical", "categorical", "categorical", "double", "categorical", "double", "categorical", "double", "categorical", "double", "categorical", "double", "double", "double", "double", "categorical", "categorical", "categorical", "categorical", "categorical", "categorical", "categorical", "double", "categorical", "categorical", "double", "double", "categorical", "categorical", "categorical", "categorical", "categorical", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
opts = setvaropts(opts, 1, "TrimNonNumeric", true);
opts = setvaropts(opts, 1, "ThousandsSeparator", ",");
opts = setvaropts(opts, [2, 4, 6, 7, 8, 10, 12, 14, 16, 21, 22, 23, 24, 25, 26, 27, 29, 30, 33, 34, 35, 36, 37], "EmptyFieldRule", "auto");
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
patn_neuropsych2 = readtable("/data/MEG/Research/ECP/Behavioural/update_pshah_060922/forR_ECP_Zscore_brief2_LisaConant.csv", opts);


%% Clear temporary variables
clear opts

%%
clear SUBNO TLESide
for i=1:size(patn_neuropsych2,1)
    SUBNO(i,:) = patn_neuropsych2.SUBNO(i);
    TLESide(i,:) = patn_neuropsych2.TLESide(i);
    Group(i,:) = patn_neuropsych2.Group(i);
end

patn_neuropsych_data = [];
patn_neuropsych_data.SUBNO = SUBNO;
patn_neuropsych_data.TLESide = TLESide;
patn_neuropsych_data.Group = Group;

