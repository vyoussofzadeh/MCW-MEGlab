% bsdir = '/data/MEG/Research/logopenicppa/BS_process';
bsdir = '/home/vyoussofzadeh/Logopenicppa_CRM';

data_tag = 'Run';
subj_ID = subj;

% tag = ['cloze',num2str(Clz)];

bsdatadir = fullfile(bsdir, 'data');
d = rdir(fullfile(bsdatadir,subj_ID,['/',data_tag,'*'],'/channel_vectorview306*.mat'));
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

%%
% switch task
%     case {1,2}
% 5: Editing channel file: run below script (while BS is open)
chan_updated = channel;
chan_updated.Channel=chan_updated.Channel(1:306);

clear a b
a = cln_data.label;
for i=1:length(chan_updated.Channel)
    b{i,:} = chan_updated.Channel(i).Name;
end
[sharedvals,idx] = intersect(b,a,'stable');
chan_updated.Channel=chan_updated.Channel(idx);
chan_updated.Projector = [];
disp('=============')
disp('channel modification was completed');

%%
savetag1 = [data_tag,'_',subj,'_',run,'_BS_channels'];
savetag2 = [data_tag,'_',subj,'_', run,'_IC_data'];

%%
Index = strfind(datafilebf{1}, '/');
tmp = datafile(Index(end)+1:end);


[data_tag1, dataf] = fileparts(datafilebf{1});

%%
% bssavedir = fullfile(bsdir,'data',subj_ID,data_tag);
bssavedir = fullfile(bsdir);
save(fullfile(bssavedir,savetag1),'chan_updated');
save(fullfile(bssavedir,savetag2),'datain','-v7.3');

%%
restoredefaultpath
addpath((allpath.ft_path));
ft_defaults
addpath(genpath(allpath.hcp_path));
addpath(genpath(allpath.cd_org));
cd(bssavedir)
