function allpath = do_init(cfg)

cd_org = cd;
addpath(genpath(cd_org));

%- FieldTrip
ft_path18 = '/MEG_data/Software/FieldTrip/fieldtrip_2018';
ft_path = fullfile(cfg.path_tools,'/fieldtrip_2022');

addpath(ft_path);
ft_defaults
ft_version

%%
allpath.path_tools = cfg.path_tools;
allpath.cd_org = cd_org;
allpath.ft_path = ft_path;
allpath.ft18 = ft_path18;

