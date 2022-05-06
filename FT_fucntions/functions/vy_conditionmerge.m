function in = vy_conditionmerge(in)

a = length(in.pst.label);
b = length(in.bsl.label);

[~,ia,ib] = intersect(in.bsl.label,in.pst.label, 'stable');

cfg = [];
if b > a
    cfg.channel = ia;
    in.bsl = ft_selectdata(cfg, in.bsl);
else
    cfg.channel = ib;
    in.pst = ft_selectdata(cfg, in.pst);
end
in.app = ft_appenddata(cfg,in.bsl,in.pst);


end
