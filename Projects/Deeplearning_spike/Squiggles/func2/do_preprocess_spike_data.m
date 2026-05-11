function [ds_data, aft] = do_preprocess_spike_data (cfg_main, datafile)

disp('preprocessing ...')
switch cfg_main.modal
    case 'meg'
        modal = 'meg';
        cfg = [];
        cfg.dataset = datafile;
        cfg.channel = {'megmag', 'meggrad', 'eog'};
        
    case 'eeg'
        modal = 'eeg';
        cfg = [];
        cfg.dataset = datafile;
        cfg.channel = {'eeg', 'eog'};
        cfg.reref = 'yes';
        cfg.refmethod = 'avg';
        cfg.refchannel = 'all';
end
raw_data = ft_preprocessing(cfg);

cfg = []; cfg.resamplefs = 500;
ds_data = ft_resampledata(cfg, raw_data);

cfg = []; cfg.toilim = [ds_data.time{:}(10*ds_data.fsample),ds_data.time{:}(end-10*ds_data.fsample)];
ds_data = ft_redefinetrial(cfg,ds_data);


% if cfg_main.applycorrection ==1
    
    %% Reject segments
    % open the browser and page through the trials
    if length(ds_data.label) > 50
        sel_sens = 1:4:length(ds_data.label);
    else
        sel_sens = 1:length(ds_data.label);
    end
    
    cfg = [];
    cfg.blocksize = ds_data.time{1}(end);
    cfg.viewmode =  'vertical'; %'butterfly';% 'vertical'; 'component'
    cfg.continuous = 'yes';
    cfg.axisfontsize = 7;
    cfg.fontsize = 7;
    cfg.preproc.demean = 'yes';
    cfg.position = [300   400   1500   400];
    cfg.channel = sel_sens;
    cfg.preproc.hpfilter = 'yes';
    cfg.preproc.hpfreq = 2;
    artf = ft_databrowser(cfg, ds_data);
    artifact_badsegment = artf.artfctdef.visual.artifact;
    
    %% Jump artifact
    cfg = [];
    cfg.artfctdef.zvalue.interactive = 'yes';
    cfg.continuous = 'yes';
    % channel selection, cutoff and padding
    cfg.artfctdef.zvalue.channel = [modal,'*'];
    
    switch cfg_main.modal
        case 'meg'
            cfg.artfctdef.zvalue.cutoff = 10;
        case 'eeg'
            cfg.artfctdef.zvalue.cutoff = 10;
    end
    cfg.artfctdef.zvalue.trlpadding = 0;
    cfg.artfctdef.zvalue.artpadding = 0;
    cfg.artfctdef.zvalue.fltpadding = 0;
    [~, artifact_jump] = ft_artifact_zvalue(cfg,ds_data);
    
    %% EOG
    cfg            = [];
    cfg.continuous = 'yes';
    % channel selection, cutoff and padding
    cfg.artfctdef.zvalue.channel     = 'EOG';
    cfg.artfctdef.zvalue.cutoff      = 6;
    cfg.artfctdef.zvalue.trlpadding  = 0;
    cfg.artfctdef.zvalue.artpadding  = 0.1;
    cfg.artfctdef.zvalue.fltpadding  = 0;
    % algorithmic parameters
    cfg.artfctdef.zvalue.bpfilter   = 'yes';
    cfg.artfctdef.zvalue.bpfilttype = 'but';
    cfg.artfctdef.zvalue.bpfreq     = [2 15];
    cfg.artfctdef.zvalue.bpfiltord  = 4;
    cfg.artfctdef.zvalue.hilbert    = 'yes';
    % feedback
    cfg.artfctdef.zvalue.interactive = 'yes';
    [~, artifact_EOG] = ft_artifact_zvalue(cfg,ds_data);
    
    %%
%     Compen_val  = 500;
%     artifact_EOG = [artifact_EOG(:,1)-Compen_val, artifact_EOG(:,1)+Compen_val];
%     artifact_jump = [artifact_jump(:,1)-Compen_val, artifact_jump(:,1)+Compen_val];
    
    aft = []; aft.eog = artifact_EOG; aft.jump = artifact_jump; aft.rejseg = artifact_badsegment;
    
%     ds_rawdata = ds_data;
    
    % if cfg_main.applycorrection ==1
    
%     fs = ds_data.fsample;
    
%     val = [];
%     val.time = ds_data.time{1};
%     val.pow = ds_data.trial{1};
%     
%     for i=1:size(aft.eog,1)
%         [~, idx1] = min(abs(val.time - aft.eog(i,1)/fs));
%         [~, idx2] = min(abs(val.time - aft.eog(i,2)/fs));
%         val.pow(:,idx1:idx2) = nan;
%     end
%     for i=1:size(aft.jump,1)
%         [~, idx1] = min(abs(val.time - aft.jump(i,1)/fs));
%         [~, idx2] = min(abs(val.time - aft.jump(i,2)/fs));
%         val.pow(:,idx1:idx2) = nan;
%     end
%     for i=1:size(aft.rejseg,1)
%         [~, idx1] = min(abs(val.time - aft.rejseg(i,1)/fs));
%         [~, idx2] = min(abs(val.time - aft.rejseg(i,2)/fs));
%         val.pow(:,idx1:idx2) = nan;
%     end
%     
%     ds_data.time{1} = val.time;
%     ds_data.trial{1} = val.pow;
% else
%     aft = [];
% end

end
%%
% cfg = [];
% cfg.artfctdef.reject = 'value'; %'partial'; %'complete'; %'value'; % this rejects complete trials, use 'partial' if you want to do partial artifact rejection
% % cfg.continuous = 'yes';
% cfg.artfctdef.eog.artifact = artifact_EOG; %
% cfg.artfctdef.jump.artifact = artifact_jump;
% cfg.artfctdef.value = 0;
% % cfg.artfctdef.muscle.artifact = artifact_muscle;
% cfg.artfctdef.xxx.artifact =
% eeg_data_cln = ft_rejectartifact(cfg,eeg_data);
% cfg.artfctdef.eog.artifact = artifact_EOG_meg; %
