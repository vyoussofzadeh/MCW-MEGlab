function [wi]  = do_time_intervals(cfg_main)

strt = cfg_main.strt; % strt = 0 sec.
spt = cfg_main.spt; % spt = 2 sec.
overlap = cfg_main.overlap; % overlap = 0.01;

wi = []; w1 = strt; l = 0.1; ov = overlap; j=1; %ov = l.*0.3
while w1+l < spt
    wi(j,:) = [w1, w1+l]; j=j+1; w1 = w1 + ov;
end
disp(wi)
% length(wi)

end