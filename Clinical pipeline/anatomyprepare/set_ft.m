restoredefaultpath
addpath(cd_org)
% ft_path = fullfile(path_tools,'/fieldtrip-20161201');
% ft_path = fullfile(path_tools,'/fieldtrip_20190419');
ft_path = '/MEG_data/Software/FieldTrip/latest/fieldtrip-master';
addpath(fullfile(path_tools,'/fieldtrip_20190419/external/spm8'));
addpath(ft_path);
ft_defaults