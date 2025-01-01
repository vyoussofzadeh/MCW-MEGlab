addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/FT_fucntions/functions_new/')


%% Fig 2 - group source source maps, whole-brain
folderPath = '/data/MEG/Research/ECP/Semantic_Decision/Results_prestim/DICS_directcontrast_localthresh/results/wb';
outputFile = 'merge_wholebrain_lateral.png';

do_combine_images(folderPath, outputFile)


%% Fig 2 - group source source maps, parcels
folderPath = '/data/MEG/Research/ECP/Semantic_Decision/Results_prestim/DICS_directcontrast_localthresh/results/parc_lat';
outputFile = 'merge_parc_lateral.png';

do_combine_images(folderPath, outputFile)

%%
folderPath = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/Main_scripts/HCP_MMP01/hcp atlas masks/HCP rois';
outputFile = 'merge_HCP_ROIs.png';

do_combine_images(folderPath, outputFile)