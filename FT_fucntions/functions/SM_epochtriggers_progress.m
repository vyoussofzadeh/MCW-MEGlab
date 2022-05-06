function [d_bsl_app, d_E_app, d_M_app, d_L_app] = SM_epochtriggers_progress(d_meg, B_time,epoch_time)

for i=1:length(B_time)
    disp([num2str(i),'/', num2str(length(B_time))])
    
    cfg = [];
    cfg.toilim = [B_time(i,1)-epoch_time, B_time(i,1)];
    data.bsl = ft_redefinetrial(cfg, d_meg);
    
    cfg = [];
    cfg.toilim = B_time(i,:);
    data.pst = ft_redefinetrial(cfg, d_meg);
    
    cfg = [];
    cfg.length = epoch_time;
    data.pst = ft_redefinetrial(cfg, data.pst);
    %%
    L = floor(length(data.pst.trial)/3);
    
    cfg = [];
    cfg.trials = 1:L;
    data.early = ft_redefinetrial(cfg, data.pst);
    cfg.trials = L+1:L*2;
    data.middle = ft_redefinetrial(cfg, data.pst);
    cfg.trials = L*2+1:L*3;
    data.late = ft_redefinetrial(cfg, data.pst);
    
    %%
    d_pst_all{i} = data.pst;
    d_bsl_all{i} = data.bsl;
    
    %-EML
    d_bsl_early{i} = data.early;
    d_bsl_middle{i} = data.middle;
    d_bsl_late{i} = data.late;
end

d_pst_app = [];
d_bsl_app = [];

for i=1:length(B_time)
    if i==1
%         d_pst_app = d_pst_all{i};
        d_bsl_app = d_bsl_all{i};
        
        d_E_app = d_bsl_early{i};
        d_M_app = d_bsl_middle{i};
        d_L_app = d_bsl_late{i};
        
    else
        cfg = [];
%         d_pst_app = ft_appenddata(cfg, d_pst_app, d_pst_all{i});
        d_bsl_app = ft_appenddata(cfg, d_bsl_app, d_bsl_all{i});
        
        %-EML
        d_E_app = ft_appenddata(cfg, d_E_app, d_bsl_early{i});
        d_M_app = ft_appenddata(cfg, d_M_app, d_bsl_middle{i});
        d_L_app = ft_appenddata(cfg, d_L_app, d_bsl_late{i});
    end
end