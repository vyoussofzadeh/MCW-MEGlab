
%% - copy data between servers
command = 'scp -r /MEG_data/MRI_database/epilepsy/DEJESUS_Elias_T1_AxBravo/nii vyoussofzadeh@squiggles.rcc.mcw.edu:/data/MEG/Vahab/';
system(command)

%% - copy data between folders
scp -r /MEG_data/Software/CAT12/cat12 /usr/local/MATLAB_Tools/spm12/toolbox
