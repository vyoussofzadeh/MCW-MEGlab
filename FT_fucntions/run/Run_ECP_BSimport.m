clear; clc, close('all'); warning off,

%%
cd '/data/MEG/Vahab/Github/MCW-MEGlab/FT';
restoredefaultpath
cd_org = cd;
addpath(genpath(cd_org));
rmpath('./Failedattemps');

indir = '/rcc/stor1/projects/ECP/MEG/MEG_Work';
%- Output dir
outdir = '/rcc/stor1/projects/ECP/MEG/';

%- Adding path
cfg_init = [];
cfg_init.path_tools = '/data/MEG/Vahab/Github/MCW-MEGlab/tools';
[allpath, atlas] = vy_init(cfg_init);

%%
load('/rcc/stor1/projects/ECP/MEG/MEG_work_BS/data/protocol.mat');

bsdir = '/rcc/stor1/projects/ECP/MEG/MEG_work_BS';
cd(bsdir)

%%
disp('1: Str')
disp('2: math');
dc = input('Select data condition:');
switch dc
    case 1
        tag = 'Str';
    case 2
        tag = 'Math';
end

%%
d = rdir([bsdir,['/*',tag,'*IC_data.mat']]);

Subj = ProtocolSubjects.Subject;
L = length(Subj);

d = rdir([indir,['/**/tSSS/*',tag,'*_raw.fif']]);
for i=1:length(d)
    [pathstr, name] = fileparts(d(i).name);
    datafolder{i} = pathstr;
    datafile{i} = d(i).name;
    Index = strfind(datafile{i}, 'EC'); Index1 = strfind(datafile{i}, '_');
    subj_all{i} = [datafile{i}(Index(2):Index1(end-3)-1)];
%     subj_all{i} = Subj(i).Name;
end

%%
set(0,'DefaultFigureWindowStyle' , 'normal')
bs_path = '/opt/matlab_toolboxes/brainstorm3';
addpath(bs_path);
brainstorm

%%
subj1 = unique(subj_all);
for ii=1:L
    [a, b] = fileparts(datafile{ii});
    if exist(fullfile(bsdir,'data',subj_all{ii},b), 'file') ~= 7
        ss = find(strcmp(subj1, subj_all{ii}), 1 );
        import_data(datafile{ii}, [], 'FT-TIMELOCK', [], ss, [], []);
        disp(ii)
    end
end

