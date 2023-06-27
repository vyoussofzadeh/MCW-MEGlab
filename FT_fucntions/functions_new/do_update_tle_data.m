function D = do_update_tle_data(cfg_main)


d_in = cfg_main.d_in;
IA = cfg_main.IA;
IB = cfg_main.IB;
var_metric = cfg_main.var_metric;
TLE_idx = cfg_main.TLE_idx;

%%
d_out = d_in;
d_out.sFiles_in = d_in.sFiles_in(IA);
d_out.sFiles_subid = d_in.sFiles_subid(IA);
d_out.all_metrics = var_metric(IB,:);

d_left = d_out;
d_left.sFiles_in = d_left.sFiles_in(TLE_idx.TLE_left);
d_left.sFiles_subid = d_left.sFiles_subid(TLE_idx.TLE_left);
d_left.all_metrics = d_left.all_metrics(TLE_idx.TLE_left,:);

d_right = d_out;
d_right.sFiles_in = d_right.sFiles_in(TLE_idx.TLE_right);
d_right.sFiles_subid = d_right.sFiles_subid(TLE_idx.TLE_right);
d_right.all_metrics = d_right.all_metrics(TLE_idx.TLE_right,:);


d_bilat = d_out;
d_bilat.sFiles_in = d_bilat.sFiles_in(TLE_idx.TLE_bilat);
d_bilat.sFiles_subid = d_bilat.sFiles_subid(TLE_idx.TLE_bilat);
d_bilat.all_metrics = d_bilat.all_metrics(TLE_idx.TLE_bilat,:);

D = [];
D.d_tle = d_out;
D.d_left = d_left;
D.d_right = d_right;
D.d_bilat = d_bilat;


