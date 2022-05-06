
if exist(savepath, 'file') == 2
    load(savepath);
else
    
    %% Artifact rejecrtion:  removing bad trials & sensors (semi-automated)
    cfg = [];
    cfg.pflag = 2; % Yes:1, No:2
    cfg.saveflag = 2; % Yes:1, No:2
    cfg.savepath = [];
    cfg.latency = [f_data.time{1}(1),f_data.time{1}(end)];
    cfg.rejectpercentage = .95;
    cfg.rbadtrl = 2;
    cfg.rbadsen = 1;
    r_data = vy_artifactreject2(cfg, f_data);
    
    %% ICA
    if flag.preprocessing.ica == 1
        cfg = [];
        cfg.savepath = outd.sub;
        cfg.savefile = savepath;
        cfg.saveflag = 1;
        cfg.overwrite = 1;
        cfg.lay = lay;
        cfg.n   = 20;
        cfg.subj = subj;
        cfg.allpath = allpath;
        cfg.select = ic_selection;
        cln_data = vy_ica_cleaning_light(cfg, r_data);
        disp('ICA cleaning was completed');
    end
end

