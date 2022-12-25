clear; clc, close('all'); warning off

%% Flags
flag.preprocessing.filtering = 1;
flag.preprocessing.artifact = 1;
flag.preprocessing.ica = 1;
flag.notch = 1;
flag.freq = 1;     % TFR & FFT
flag.time = 1;     % Time-locked & Cov estimation
flag.gave = 0;     % grand average analysis
flag.anatomy = 1;     % grand average analysis
flag.sourceanalysis = 1;     % grand average analysis
flag.speechanalysis = 2;     % speech analysis

%% Initial settings
% set(0,'DefaultFigureWindowStyle','normal')

cd '/MEG_data/Vahab/Github/MCW-MEGlab/FT';
restoredefaultpath
cd_org = cd;
addpath(genpath(cd_org));

%- Input dir
indir = '/MEG_data/epilepsy';
%- Output dir
outdir = '/MEG_data/Vahab/Processed_data';

%- Adding path
cfg_init = [];
cfg_init.path_tools = '/MEG_data/Vahab/Github/tools';
[allpath, atlas] = vy_init(cfg_init);

%%
disp('1: Definition naming')
disp('2: Picture naming');
disp('3: Motor');
disp('4: Somatosensory');
task = input('Eneter the task: ');
switch task
    case 1
        %- Auditory definition naming
        tag = 'DFN';
    case 2
        %- Visual picture naming
        tag = 'PN';
    case 3
        %- Motor task
        tag = 'motor';
    case 4
        %- somatosensory median nerve stim
        tag = 'SSEF';
end

%%
cd(indir)
[subjdir] = uigetdir;

%%
d = rdir([subjdir,['/**/','sss','/*',tag,'*/*raw_tsss.fif']]);

switch task
    case 4
        d = rdir([subjdir,['/**/','sss','/*',tag,'*/*tsss.fif']]);
end

%%
clear subj datafolder datafile datafile1
for i=1:length(d)
    [pathstr, name] = fileparts(d(i).name);
    datafolder{i} = pathstr;
    datafile{i} = d(i).name;
    Index = strfind(datafile{i}, '/');
    subj = datafile{i}(Index(3)+1:Index(4)-1);
end
datafile1 = datafile';
disp(datafile1)
if length(datafile1) > 1
    datasel = input('choose data to analyze, eg, 1,2:');
else
    datasel = 1;
end
disp([subj, ' and,'])
disp(datafile1{datasel})
disp('was selected for the analysis.')
disp('============');

%%
% epoch_type = 'STI101';
% 
% %% 4D layout
% cfg = [];
% cfg.layout = 'neuromag306mag.lay';
% lay = ft_prepare_layout(cfg);
% % ft_layoutplot(cfg);
% disp('============');
% 
% %%
% close all
% datafile = datafile1{datasel}; % spm_select(inf,'dir','Select MEG folder'); % e.g. H:\VNS\MEG\C-105\CRM\1
% Index = strfind(datafile, '/');
% Date  = datafile(Index(4)+1:Index(5)-1);
% disp('============');
% disp(datafile)
% disp(['subj:',subj])
% disp(['Date:',Date])
% disp('============');
% 
% %%
% %-elec/grad
% sens = ft_read_sens(datafile);
% sens = ft_convert_units(sens,'mm');
% 
% %%
% outd.sub = fullfile(outdir,'ft_process',subj, tag);
% if exist(outd.sub, 'file') == 0
%     mkdir(outd.sub);   %create the directory
% end
% cd(outd.sub)
% disp(['outputdir:',outd.sub])
% disp('============');

%%
bsdatadir = fullfile(indir,subj,'brainstorm_db/data');
d = rdir(fullfile(bsdatadir,subj,['/@*',tag,'*'],'/channel_vectorview306*.mat'));
if isempty(d)
    d = rdir(fullfile(bsdatadir,subj,['/*',tag,'*'],'/channel_vectorview306*.mat'));
    if exist('tag1','var') && isempty(d)
        d = rdir(fullfile(bsdatadir,subj,['/*',tag1,'*'],'/channel_vectorview306*.mat'));
    end
end

%%
for i=1:length(d)
    [pathstr, name] = fileparts(d(i).name);
    %     databf{i} = pathstr;
    datafilebf{i} = d(i).name;
end
disp(datafilebf')
channel = load(datafilebf{1});
disp('=============')
disp(datafilebf{1})
disp('was loaded')

% round(1e3.*channel.HeadPoints.Loc(:,1:3))'
% channel.HeadPoints.Label(1:3)

disp([channel.HeadPoints.Label(1)])
disp(1e3.*channel.HeadPoints.Loc(:,1)')
disp([channel.HeadPoints.Label(2)])
disp(1e3.*channel.HeadPoints.Loc(:,2)')
disp([channel.HeadPoints.Label(3)])
disp(1e3.*channel.HeadPoints.Loc(:,3)')

% T = [-0.0278770366090040,0.999611359894383,-3.93522432922975e-08,7.35523790542143e-05;-0.999611359894383,-0.0278770366090059,-4.90540103330250e-08,0.00263796671944081;-5.01319699041742e-08,3.79694689904360e-08,0.999999999999998,3.85758759629082e-09;0,0,0,1];
% 
% hc_tr = ft_transform_geometry(T,channel.HeadPoints.Loc);
sMRI = load(sMRI1);
P_mri = cs_convert(sMRI, 'scs', 'mri',(1e3.*channel.HeadPoints.Loc(:,2))')

%% Preprocesssing
% clear cln_data
% ic_selection = 1; % 1: manual, 2: automated
% 
% hdr = ft_read_header(datafile);
% 
% addpath('/MEG_data/Vahab/Github/MCW-MEGlab/FT/Clinical pipeline');
% [hc] = read_neuromag_hc(datafile);
% 
% % hdr = ft_read_header(datafile, 'checkmaxfilter', false);
% % nFid = size(hdr.orig.dig,2);
% %%
% clc
% disp(hc.dewar.label(1:3));
% disp(hc.dewar.pos(1:3,:));
% 
% % fid_neuromag = hc.dewar.pos(1:3,:);
% 
% %%
% % From Chan files
% transform = [-0.0278770366090040	0.999611359894383	-3.93522432922975e-08	7.35523790542143e-05
% -0.999611359894383	-0.0278770366090059	-4.90540103330250e-08	0.00263796671944081
% -5.01319699041742e-08	3.79694689904360e-08	0.999999999999998	3.85758759629082e-09
% 0	0	0	1];
% 
% hc_tr = ft_transform_geometry(transform, hc.dewar);
% % clc
% disp(hc_tr.label(1:3)');
% disp(hc_tr.pos(1:3,:));
% 
% %%
mridir = fullfile(indir,subj,'brainstorm_db/anat');
d = rdir(fullfile(mridir,subj,'subjectimage*.mat'));
sMRI1 = d.name;
load(sMRI1);
% fid.SCS = SCS;
% fid.NCS = NCS;
% A = load(sMRI1);

brainstorm
P_mri = cs_convert(sMRI1, 'voxel', 'mri', channel.HeadPoints.Loc(:,1));   % Voxel => MRI coordinates


