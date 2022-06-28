

ChannelMat = in_bst_channel('EC1002/ec1002_SD_run2_raw_clean/channel_vectorview306_acc1.mat');
HeadModelFile = 'EC1002/ec1002_SD_run2_raw_clean/headmodel_surf_os_meg.mat';
HeadModelMat = in_bst_headmodel(HeadModelFile);
[ftHeadmodel, ftLeadfield, iChannelsData] = out_fieldtrip_headmodel(HeadModelMat, ChannelMat, 1:306, 1);

grad = [];


dataft = load('/data/MEG/Research/spendl/Shared_scripts/datafiles/test_ft.mat')


% Interpolate bad channels: Find neighbours
cfg             = [];
cfg.method      = 'triangulation';
cfg.senstype   = 'MEG'; %
cfg.grad       = dataft.grad;
neighbours_MEG = ft_prepare_neighbours(cfg, dataft);


% plotting neighbours for inspection
cfg            = [];
cfg.neighbours = neighbours_MEG;
cfg.senstype   = 'MEG';
ft_neighbourplot(cfg, dataft);

% Interpolate channels
cfg = [];
cfg.method                    = 'spline';
cfg.neighbours                = neighbours_MEG;
cfg.badchannel                = badchannels;
cfg.senstype                  = 'MEG';
interpolated_data = ft_channelrepair(cfg, dataft);