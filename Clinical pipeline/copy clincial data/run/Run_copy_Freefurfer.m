
%% Find FS files
fs_files = '/MEG_data/LAB_MEMBERS/Vahab/Github/MCW-MEGlab/FT/Clinical pipeline/copy clincial data/fs_details.mat';
if ~~exist(fs_files,'file')
    load(fs_files)
else
    %- look up the FS from brainstorm history
    clc
    idx = [];
    idx.matched = data_idx(ib);
    
    clc
    clear datafile_sub
    clear fs_hist
    kk=1; km=1;
    
    for i=1:length(name_list.newmatched)
        fname = name_list.newmatched{i};
        idx.matched2 = find(contains(name_list.all,fname)==1);
        
        if ~isempty(idx.matched2)
            savepfolder = fullfile(savedir, name_list.all{idx.matched2}, 'anat');
            if exist(savepfolder, 'file') == 0, mkdir(savepfolder), end
            cd(savepfolder)
            d = dir(fullfile(indir,fname, '/brainstorm_db'));
            if ~isempty(d)
                cd(fullfile(indir,fname))
                d_anat = rdir(fullfile(indir,fname,'/**/', fname, 'subjectimage_*.mat'));
                d = d_anat;
                
                if ~isempty(d)
                    if length(d) > 1
                        disp('------')
                        for ii=1:length(d)
                            disp([num2str(ii), ': ', d(ii).name])
                        end
                        sel = input('sel file:');
                    else
                        sel = 1;
                    end
                    [pathstr, name] = fileparts(d(sel).name);
                    a = load(d(sel).name);
                    if strcmp(a.Comment,'MRI') && isfield(a, 'History')
                        fs_hist{kk} = a.History{1,3}(14:end);
                        disp(fs_hist{kk}); kk=kk+1;
                    end
                else
                    disp(['missing subjectimage_T1.mat for ', fname])
                    missing_data(km) = i; km = km+1;
                end
            end
        end
    end
    
    % Get FS subject details
    clc
    fs_hist_sub = [];
    for i=1:length(fs_hist)
        disp([num2str(i),'/',num2str(length(fs_hist))]);
        
        tkz = tokenize(fs_hist{i},'/');
        disp(tkz{5})
        sub = find(contains (name_list.all,tkz{5},'Ignorecase', true) ==1);
        
        if ~isempty(sub) && length(sub) ==1
            disp(name_list.all{sub})
            sub_sel = 1;
        else
            disp('no match was found')
            sub_sel = 2;
        end
        disp('1:yes')
        disp('2:no')
        
        switch sub_sel
            case 1
                ss = name_list.all{sub};
            case 2
                disp(name_list.all)
                tkz = tokenize(fs_hist{i},'/');
                if contains(fs_hist{i},'/epilepsy')
                    idx = strfind(fs_hist{i},'/epilepsy');
                    idx2 = strfind(fs_hist{i}(idx+10:end),'/');
                    string = fs_hist{i}(idx+10:idx+10+idx2-2);
                else
                    string = tkz{end-2};
                end
                %- IN PROGRESSl
                TF = find(startsWith(name_list.all,lower(string(1)))==1);
                name_list.all(TF)
                disp(fs_hist{i});
                ss = input('enter sub name:','s');
        end
        fs_hist_sub{i} = ss;
    end
    save(fs_file,'fs_hist_sub','fs_hist')
end

%%
for i=1:length(fs_hist_sub)
    savefolder = fullfile(savedir,fs_hist_sub{i},'/anat/');
    if exist(savefolder, 'dir') == 0
        mkdir(savefolder);   %create the directory
    end
    cd(savefolder)
    dcheck = rdir([savefolder,'/*/*']);
    dcheck2 = rdir([savefolder,'/*.tar']);
    if exist(fs_hist{i},'file') && isempty(dcheck) && isempty(dcheck2)
        Index = strfind(fs_hist{i}, '/');
        if size(Index,2) >= 6
            disp(fs_hist{i}(1:Index(6)-1)); disp('to'),
            disp(savefolder)
            copyfile(fs_hist{i}(1:Index(6)-1),savefolder),
        else
            disp(fs_hist{i})
            fs_innew = input('enter better FS folder dir:');
            copyfile(fs_innew,savefolder),
        end
    elseif ~exist(fs_hist{i},'file') && isempty(dcheck)
        disp(fs_hist{i})
        disp('does not exist!')
        d_tar = rdir([fs_hist{i}(1:end-6),'*.tar']);
        dcheck = rdir([savefolder,'/*']);
        if ~isempty(d_tar)  && isempty(dcheck)
            disp(d_tar.name)
            disp(d_tar.name); disp('to'),
            disp(savefolder)
            copyfile(d_tar.name,savefolder),
        elseif  isempty(d_tar)  && isempty(dcheck)
            fs_innew = input('enter better FS folder dir:','s');
            disp(fs_innew); disp('to'),
            disp(savefolder)
            copyfile(fs_innew,savefolder),
        end
        disp('data copied')
    end
end

%% Copy the remiander FS files, manually
missing_sub = [];
for i=1:length(missing_data)
    missing_sub{i} = name_list.newmatched{missing_data(i)};
end

clc
for i=1:length(missing_sub)
    disp([num2str(i),'/',num2str(length(missing_sub))]);
    
    savefolder = fullfile(savedir,missing_sub{i},'/anat/');
    if exist(savefolder, 'dir') == 0, mkdir(savefolder);   end
    cd(savefolder)
    dcheck = rdir([savefolder,'/*']); dcheck1 = rdir([savefolder,'*/*']);
    if isempty(dcheck) && isempty(dcheck1)
        disp(missing_sub{i})
        
        % FS file
        tkz = tokenize(missing_sub{i},'_');
        d_fs = dir(fullfile('/MEG_data/MRI_database/epilepsy/',...
            [upper(tkz{1}), '_', upper(tkz{2}(1)),tkz{2}(2:end)],'/*FSRecon*'));
        
        % tar file
        d_tar = rdir(['/MEG_data/MRI_database/epilepsy/',...
            [upper(tkz{1}), '_', upper(tkz{2}(1)),tkz{2}(2:end)],'*/*.tar']);
        
        if ~isempty(d_tar) && size(d_tar,1) ==1
            copyfile(d_tar.name,savefolder),
        elseif ~isempty(d_fs) && size(d_fs,1) ==1
            copyfile(fullfile(d_fs.folder, d_fs.name),savefolder)
        else
            fs_innew = input('enter better FS folder dir:','s');
            copyfile(fs_innew,savefolder),
        end
        disp('data was copied')
    end
end

