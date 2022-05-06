function [data_pst_app,data_bsl_app] = SM_epochtriggers(data_meg, B_time,epoch_time)

for i=1:length(B_time)
    disp([num2str(i),'/', num2str(length(B_time))])
    
    cfg = [];
    cfg.toilim = [B_time(i,1)-epoch_time, B_time(i,1)];
    data.bsl = ft_redefinetrial(cfg, data_meg);
    
    cfg = [];
    cfg.toilim = B_time(i,:);
    data.pst = ft_redefinetrial(cfg, data_meg);
    
    cfg = [];
    cfg.length = epoch_time;
    data.pst = ft_redefinetrial(cfg, data.pst);
    
    %- replicating the baseline (to have identical trials as post-stim)
%     tmp = data.pst;
%     for j=1:length(data.pst.trial)
%         tmp.trial{j} = data.bsl.trial{1};
%         tmp.time{j} = data.bsl.time{1};
%     end
%     data.bsl = tmp;
    
    data_pst_all{i} = data.pst;
    data_bsl_all{i} = data.bsl;
end

data_pst_app = [];
data_bsl_app = [];

for i=1:length(B_time)
    if i==1
        data_pst_app = data_pst_all{i};
        data_bsl_app = data_bsl_all{i};
    else
        cfg = [];
        data_pst_app = ft_appenddata(cfg, data_pst_app, data_pst_all{i});
        data_bsl_app = ft_appenddata(cfg, data_bsl_app, data_bsl_all{i});
    end
end