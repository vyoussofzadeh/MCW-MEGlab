savepath = fullfile(outd.sub,['trial_info_',subj,'.mat']);
if exist(savepath, 'file') == 2
    load(savepath);
else
    %%
    trlInfoColDescr = [];
    
    scanmnem = 'StoryM_neuromag';
    trialFunName=['trialfun_',scanmnem];
    trialDefFuncFile=['alltrialdefparams_',scanmnem];
    trialDefFuncHandle=str2func(trialDefFuncFile);
    savesuffix_trialinfo='trialinfo';
    
    experimentid = 'SM';
    scanid = run;
    pipelineprefix='baddata';
    
    resultprefix = sprintf('%s_%s_%s', experimentid, scanid, pipelineprefix);
    %------------------------------------------
    %--- Extract trials definition
    allTrlCfgs=trialDefFuncHandle();
    Ncasegroups=length(allTrlCfgs);
    outputTrialInfoSummaryFile = [resultprefix,'_raw',savesuffix_trialinfo,'_QC.txt'];
    montage         = hcp_exgmontage(subj, experimentid, scanid);
    iCaseGroup=2;
%     iCaseGroup=2;
    
    %------------------------------------------
    trlCfg                 = allTrlCfgs{iCaseGroup};
    trlCfg.datafile        = datafile;
    trlCfg.trialfun        = trialFunName;
    if ~isempty(regexp(scanid,'Motor', 'once'))
        trlCfg.trialdef.montage=montage;
    end
    %------------------------------------------
    eval(['[trl,trialinfo,trlInfoColDescr,trialSummary,scanStartSamp,scanEndSamp,warninfo]=',trialFunName,'(trlCfg);']);
    % Save the summary ascii file
    hcp_write_ascii(outputTrialInfoSummaryFile,'trialSummary'); % Trial Summary should be the same no matter what the input arguments are. This is because the trial definition function creates all information about the trial definition. This is what the summary contains. The trl field only contains the trials defined by the input arguments.
    ErrorFile=['ERROR_',outputTrialInfoSummaryFile];
    delete([ErrorFile,'*']); % Remove any error files present from previous runs
    WarningFile=['WARNING_',outputTrialInfoSummaryFile];
    if isempty(warninfo)
        delete([WarningFile,'*']); % Remove any warning files present from previous runs
    else
        hcp_write_ascii(WarningFile,'warninfo');
    end
    save(savepath, 'trl', 'trlInfoColDescr','trialSummary','trialinfo');
end