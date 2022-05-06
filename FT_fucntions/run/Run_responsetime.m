savefile = [tsk,'_RT.mat'];
if exist(savefile, 'file') == 2
    load(savefile)
else
    clear RT
    tt_all = []; k = 1;
    for i=1:length(files_sel)
        load(files_sel{i});
        datafile = files_sel{i}; % spm_select(inf,'dir','Select MEG folder'); % e.g. H:\VNS\MEG\C-105\CRM\1
        Index = strfind(datafile, '/');
        Date = datafile(Index(5)+1:Index(6)-1);
        Subj  = datafile(Index(6)+1:Index(7)-1);
        RT(i,1) = mtt; RT(i,2) = stt;
        disp(datafile)
        disp(['subj:',Subj])
        disp(['Date:',Date])
        disp([num2str(i),'/', num2str(length(files_sel))])
        disp('============');
        names{i} = [num2str(i),'-',Subj,'-',Date];
        Sub_all{i} = Subj;
        Date_all{i}= Date;
    end
    RT_all = table(Sub_all',Date_all', RT);
    disp(RT_all)
    
    figure,plot(RT(:,1),'*')
    save([tsk,'_RT.mat'],'RT_all', 'RT','-v7.3');
end



