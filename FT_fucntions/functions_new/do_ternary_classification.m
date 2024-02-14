function ternary_values = do_ternary_classification(cfg)

LI = cfg.LI;
thre = cfg.thre;


% LI = (LI - min(LI(:))) ./ (max(LI(:)) - min(LI(:)));

LI = LI./max(LI(:));


ternary_values = zeros(size(LI)); % Preallocate for speed

for i=1:size(LI,1)
    for j = 1:size(LI,2)
        % Generate ternary values
        ternary_values(i,j, LI(i,j,:) < -1*thre) = -1;
        ternary_values(i,j,LI(i,j,:) > thre) = 1;
    end
end
end
