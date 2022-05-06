
tag = ['cloze',num2str(Clz)];

datacheck = fullfile(bsdir,'data',subj,[tag,'_',subj,'_', stag, '*'],'*.mat');
d = rdir(datacheck);

savetag2 = [tag,'_',subj,'_', stag,'_IC_data.mat'];
% if exist(fullfile(bsdir,savetag2), 'file') ~= 2
if isempty(d)
    
    bsdatadir = fullfile(bsdir, 'data');
    d = rdir(fullfile(bsdatadir,subj,['/',tag,'*'],'/channel_vectorview306*.mat'));
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
    if ~isempty(d)
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
        savetag1 = [tag,'_',subj,'_',stag,'_BS_channels'];
        savetag2 = [tag,'_',subj,'_', stag,'_IC_data'];
        
        %%
        bssavedir = fullfile(bsdir,'ftExport');
        if exist(bssavedir, 'file') == 0, mkdir(bssavedir);   end
        
        save(fullfile(bssavedir,savetag1),'chan_updated');
        save(fullfile(bssavedir,savetag2),'datain','-v7.3');
        
        %%
        restoredefaultpath
        addpath((allpath.ft_path));
        ft_defaults
        addpath(genpath(allpath.hcp_path));
        addpath(genpath(allpath.cd_org));
        cd(bssavedir)
    end
end