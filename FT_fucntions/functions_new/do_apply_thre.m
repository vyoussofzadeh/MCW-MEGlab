function Val_thre = do_apply_thre(cfg, In_val)

if cfg.do_atan ==1
    In_val = atan(In_val);
end

if cfg.do_normal ==1
    In_val = (In_val - min(In_val(:))) ./ (max(In_val(:)) - min(In_val(:)));
end

idx = In_val >= cfg.thre.*max(In_val([cfg.idx_L, cfg.idx_R]));
Val_thre = zeros(size(In_val)); Val_thre(idx) = In_val(idx);


if cfg.do_plot ==1
    figure, plot(In_val), hold on, yline(cfg.thre), plot(Val_thre,'Color',[0.5 0.5 0.5])
end