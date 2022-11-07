%- Copying MEG files to a temporary folder
clc
idx = [];
idx.matched = data_idx(ib);

clc
clear datafile_sub
k=1;
for i=1:length(name_list.newmatched)
    
    fname = name_list.newmatched{i};
    idx.matched2 = find(contains(name_list.all,fname)==1);
    
    if ~isempty(idx.matched2)
        
        savepfolder = fullfile(savedir, name_list.all{idx.matched2}, 'sss');
        
        d1 = rdir(fullfile(indir,fname, '/**/sss/**/*spont*tsss*raw.fif'));
        d2 = rdir(fullfile(indir,fname, '/**/sss/**/*spont*sss_ecgClean*raw.fif'));
        d3 = rdir(fullfile(indir,fname, '/**/sss/**/*spont*t_sss*raw.fif'));
        d4 = rdir(fullfile(indir,fname, '/**/sss/**/*emptyroom*sss*raw.fif'));
        d5 = rdir(fullfile(indir,fname, '/**/sss/**/*spont*tsss*raw.fif'));
        d6 = rdir(fullfile(indir,fname, '/**/sss/**/*spont*sss*raw.fif'));
        d7 = rdir(fullfile(indir,fname, '/**/sss/**/*spont*raw_tsss.fif'));
        d_all = [d1; d2; d3; d4; d5; d6; d7];
        
        cd(fullfile(indir,fname))
        
        d = d_all;
        clear subj_all datafile
        for ii = 1:length(d)
            [pathstr, name] = fileparts(d(ii).name);
            datafolder{ii} = pathstr;
            datafile{ii} = d(ii).name;
            Index = strfind(datafile{ii}, '/');
            subj_all{ii} = datafile{ii}(Index(5)+1:Index(6)-1);
            %         disp(datafile{ii})
        end
        datafile_sub{k} = datafile';
        k=1+k;
        
        if exist(savepfolder, 'file') == 0, mkdir(savepfolder), end
        
        for j=1:length(datafile')
            tkz = tokenize(datafile{j},'/');
            savepfile = fullfile(savepfolder, tkz{end});
            if ~exist(savepfile,'file')
                disp([num2str(i),': copying: ', fname])
                disp(datafile{j})
                disp('to')
                disp(savepfile)
                %                 pause
                copyfile(datafile{j}, savepfile),
            end
        end
    else
        disp([num2str(i), ': no detected: ', fname])
    end
end