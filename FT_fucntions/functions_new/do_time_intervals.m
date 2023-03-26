function [wi]  = do_time_intervals(cfg)

wi = []; w1 = 0; l = 0.1; ov = 0.01; j=1; %ov = l.*0.3
while w1+l < 2
    wi(j,:) = [w1, w1+l]; j=j+1; w1 = w1 + ov;
end
disp(wi)
% length(wi)

end