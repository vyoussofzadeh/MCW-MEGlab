function cln_data = vy_ica_cleaning_light_SZ(cfg_main, data_in)

% satis = 0;
disp('ica cleaning ...');
if exist(cfg_main.savefile, 'file') && (cfg_main.overwrite) == 2
    load(cfg_main.savefile)
else
    
    savepath = fullfile(cfg_main.savepath,'ica');
    if exist(savepath, 'file') == 0, mkdir(savepath), end
    
    n = cfg_main.n; % ICs
    
    cfg = [];
    cfg.lay = cfg_main.lay;
    cfg.subj = cfg_main.subj;
    cfg.n = n;
    cfg.savefig = 1;
    cfg.allpath = cfg_main.allpath;
    comp = vy_ica(cfg, data_in);
    %     title(savepath)
    cfg = [];
    cfg.updatesens = 'no';
    
    
    %% Rejecting bad ICAs and prjecting back into data space
    disp('=============================')
    cfg = [];
    if cfg_main.select == 1
        bic = input(['Select bad ICs for ' cfg_main.subj,':']);
        cfg.component = comp.label(bic);
    else
        cfg.component =  report.bchan;
    end
    cfg.updatesens = 'no';
    cln_data = ft_rejectcomponent(cfg, comp, data_in);
    close all
    if cfg_main.saveflag == 1
        save(cfg_main.savefile, 'cln_data', '-v7.3');
        textfile_rej = 'ica/selected_badICs';
        badICs = cell2table(cfg.component);
        if ~isempty(badICs)
            badICs.Properties.VariableNames{'Var1'} = 'bICAs';
            writetable(badICs,textfile_rej,'Delimiter',' ');
        end
    end
    
    %% back to new ft
    restoredefaultpath
    addpath((cfg_main.allpath.ft_path));
    ft_defaults
    addpath(genpath(cfg_main.allpath.hcp_path));
    addpath(genpath(cfg_main.allpath.cd_org));
    addpath(genpath(cfg_main.allpath.exfig_path));
    
end

