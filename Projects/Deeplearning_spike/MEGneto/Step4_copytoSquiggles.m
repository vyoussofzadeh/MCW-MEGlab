%l% Copy to Squiggles
cd ('/MEG_data/Research_studies/Epil_annotated_data')
command = 'scp -r /MEG_data/Research_studies/Epil_annotated_data/annotated_data vyoussofzadeh@squiggles.rcc.mcw.edu:/data/MEG/Research/SpikeDectection/Epil_annotated_data';
system(command)


command = 'scp -r /MEG_data/Research_studies/Epil_annotated_data/annotated_data_nospike vyoussofzadeh@squiggles.rcc.mcw.edu:/data/MEG/Research/SpikeDectection/Epil_annotated_data';
system(command)


command = 'scp -r /MEG_data/Research_studies/Epil_annotated_data/annotated_info vyoussofzadeh@squiggles.rcc.mcw.edu:/data/MEG/Research/SpikeDectection/Epil_annotated_data';
system(command)