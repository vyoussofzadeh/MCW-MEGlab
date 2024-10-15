addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/FT_fucntions/functions_new/')


folderPath = '/data/MEG/Research/ECP/Semantic_Decision/Results_prestim/DICS_directcontrast_localthresh/results/wb';
outputFile = 'merge_wholebrain_lateral.png'

do_combine_images(folderPath, outputFile)


folderPath = '/data/MEG/Research/ECP/Semantic_Decision/Results_prestim/DICS_directcontrast_localthresh/results/parc_lat';
outputFile = 'merge_parc_lateral.png'

do_combine_images(folderPath, outputFile)

