function data_app = vy_merge_data(data1, data2)

a = length(data1.label); b = length(data2.label);

% tmp1 = data1;
% tmp2 = data2;

[~,ia,ib] = intersect(data2.label,data1.label, 'stable');

cfg = [];
if b > a
    cfg.channel = ia; data2 = ft_selectdata(cfg, data2);
    cfg.channel = ib; data1 = ft_selectdata(cfg, data1);
else
    cfg.channel = ib; data2 = ft_selectdata(cfg, data2);
    cfg.channel = ia; data1 = ft_selectdata(cfg, data1);
end

cfg.keepsampleinfo = 'yes';
data_app = ft_appenddata(cfg, data2, data1);
data_app.hdr = data2.hdr;