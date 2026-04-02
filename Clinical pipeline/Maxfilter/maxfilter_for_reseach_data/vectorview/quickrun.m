

%% terminal command
mkdir -p /tmp/maxfilter_test && /neuro/bin/util/maxfilter -f "/MEG_acq/cnrp_tacs_healthy/tsss_move_done/hc003_v2/220810/orig/PSTM_Block1_raw.fif" -o "/tmp/maxfilter_test/PSTM_Block1_raw.fif" -ctc /neuro/databases/Elekta_Orig_Files/ctc/ct_sparse.fif -cal /neuro/databases/Elekta_Orig_Files/sss-oldElekta/sss_cal.dat -autobad off -st 10 -corr .9 -trans "/MEG_acq/cnrp_tacs_healthy/tsss_move_done/hc003_v2/220810/orig/PSTM_Block1_raw.fif" -force