


alpha  = 11;                       % chose new p-cut-off
mask   = (test.tmap > alpha);        % significant samples
t_thr  = test.tmap .* mask;          % masked t-values

test2 = test;
test2.tmap = t_thr;