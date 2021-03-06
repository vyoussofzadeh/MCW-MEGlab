% Script generated by Brainstorm (19-Aug-2019)

% Input files
sFiles = {...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_02.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_03.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_04.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_05.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_06.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_07.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_08.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_09.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_10.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_100.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_101.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_102.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_11.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_12.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_13.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_14.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_15.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_16.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_17.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_18.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_19.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_20.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_21.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_22.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_23.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_24.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_25.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_26.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_27.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_28.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_29.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_30.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_31.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_32.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_33.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_34.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_35.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_36.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_37.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_38.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_39.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_40.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_41.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_42.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_43.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_44.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_45.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_46.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_47.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_48.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_49.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_50.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_51.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_52.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_53.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_54.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_55.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_56.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_57.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_58.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_59.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_60.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_61.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_62.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_63.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_64.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_65.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_66.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_67.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_68.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_69.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_70.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_71.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_72.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_73.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_74.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_75.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_76.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_77.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_78.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_79.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_80.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_81.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_82.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_83.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_84.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_85.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_86.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_87.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_88.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_89.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_90.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_91.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_92.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_93.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_94.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_95.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_96.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_97.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_98.mat', ...
    'bednar_peggy/DFN_bednar_peggy_IC_data/data_DFN_bednar_peggy_IC_99.mat'};

% Start a new report
bst_report('Start', sFiles);

% Process: FieldTrip: ft_sourceanalysis
sFiles = bst_process('CallProcess', 'process_ft_sourceanalysis_vahab', sFiles, [], ...
    'method',     'dics', ...  % DICS beamformer
    'sensortype', 'MEG');  % MEG

% Save and display report
ReportFile = bst_report('Save', sFiles);
bst_report('Open', ReportFile);
% bst_report('Export', ReportFile, ExportDir);

