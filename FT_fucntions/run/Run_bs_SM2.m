% bsdir = '/rcc/stor1/projects/ECP/MEG/MEG_work_BS';
% tag = ['cloze',num2str(Clz)];
cd(bsdir)
savetag2 = [tag.dcon,'_',subj,'_', tag.task,'_', tag.dg, '_IC_data.mat'];
if exist(fullfile(bsdir,savetag2), 'file') ~= 2
    
    %%
    bsdatadir = fullfile(bsdir, 'data');
    d = rdir(fullfile(bsdatadir,subj,['/*',tag.task,'*'],'/channel_vectorview306*.mat'));
    if isempty(d)
        d = rdir(fullfile(bsdatadir,subj,['/*',tag.dcon,'_', tag.dg,'*'],'/channel_vectorview306*.mat'));
        if exist('tag1','var') && isempty(d)
            d = rdir(fullfile(bsdatadir,subj,['/*',tag1,'*'],'/channel_vectorview306*.mat'));
        end
    end
    
    %%
    clear nn
    for i=1:length(d)
        [pathstr, name] = fileparts(d(i).name);
        %     databf{i} = pathstr;
        datafilebf{i} = d(i).name;
        nn{i} = [num2str(i),':',d(i).name];
    end
    % disp(datafilebf')
    if ~isempty(d)
        
        disp(nn')
        if length(nn)>1
%             sel = input('choose data?');
%             sel = length(nn);
            for dpick=1:length(nn)
                Index = strfind(datafilebf{dpick}, 'raw');
                if ~isempty(Index), sel = dpick; break, end
            end
        else
            sel = 1;
        end
        channel = load(datafilebf{sel});
        disp('=============')
        disp(datafilebf{sel})
        disp('was loaded')
        
        %%
        % switch task
        %     case {1,2}
        % 5: Editing channel file: run below script (while BS is open)
        chan_updated = channel;
        chan_updated.Channel = chan_updated.Channel(1:306);
        
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
        savetag1 = [tag.dcon,'_',subj,'_',tag.task,'_', tag.dg, '_BS_channels'];
        savetag2 = [tag.dcon,'_',subj,'_', tag.task,'_', tag.dg, '_IC_data'];
        
        %%
        % bssavedir = fullfile(bsdir,'data',subj);
        save(fullfile(bsdir,savetag1),'chan_updated');
        save(fullfile(bsdir,savetag2),'datain','-v7.3');
        
        %%
        restoredefaultpath
        addpath((allpath.ft_path));
        ft_defaults
        addpath(genpath(allpath.hcp_path));
        addpath(genpath(allpath.cd_org));
        cd(bsdir)
        
    end
end