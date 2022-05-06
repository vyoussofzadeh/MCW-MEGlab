%%
epoch_time = 2;
switch dc
    case 1
        %- Story
        B_str_time = (B_str/sample_scale)/data_raw_resampled.fsample;
        [data_str_pst,data_str_bsl] = SM_epochtriggers(data_raw_resampled, B_str_time,epoch_time);
        
    case 2
        %- Math
        B_math_time = (B_math/sample_scale)/data_raw_resampled.fsample;
        [data_math_pst,data_math_bsl] = SM_epochtriggers(data_raw_resampled, B_math_time, epoch_time);
    case 3
        %- Story-Math
        B_str_time = (B_str/sample_scale)/data_raw_resampled.fsample;
        [data_str_pst,data_str_bsl] = SM_epochtriggers(data_raw_resampled, B_str_time,epoch_time);
        B_math_time = (B_math/sample_scale)/data_raw_resampled.fsample;
        [data_math_pst,data_math_bsl] = SM_epochtriggers(data_raw_resampled, B_math_time, epoch_time);
end

%%
switch dc
    case 1
        cln_str_pst = SM_preprocess(data_str_pst,lay,subj, allpath);
        cln_str_bsl = SM_preprocess(data_str_bsl,lay,subj, allpath);
    case 2
        cln_math_pst = SM_preprocess(data_math_pst,lay,subj, allpath);
        cln_math_bsl = SM_preprocess(data_math_bsl,lay,subj, allpath);
    case 3
        cln_str_pst = SM_preprocess(data_str_pst,lay,subj, allpath);
        cln_math_pst = SM_preprocess(data_math_pst,lay,subj, allpath);
end

%%
switch dc
    case 1
        p_data_str_math.pst = SM_adjust_time(cln_str_pst,epoch_time);
        p_data_str_math.bsl = SM_adjust_time(cln_str_bsl,epoch_time);
    case 2
        p_data_str_math.pst = SM_adjust_time(cln_math_pst,epoch_time);
        p_data_str_math.bsl = SM_adjust_time(cln_math_bsl,epoch_time);
    case 3
        p_data_str_math.pst = SM_adjust_time(cln_str_pst,epoch_time);
        p_data_str_math.bsl = SM_adjust_time(cln_math_pst,epoch_time);
end