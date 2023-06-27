function [IA, IB, var_metric, var_name, TLE_idx] = do_extract_tle_data(cfg_main)


patn_neuropsych_tle = cfg_main.patn_neuropsych_tle;
TLESide = patn_neuropsych_tle.TLESide; 
SUBNO = patn_neuropsych_tle.SUBNO;
AEDcount = patn_neuropsych_tle.AEDcount;
WASI_BlckR_ZScore = patn_neuropsych_tle.WASI_BlckR_ZScore;
sub_pt = cfg_main.sub_pt;

%%
k = 1;
var_name = [];
var_metric = [];
for i= [47:94] %[10, 14, 16, 18, 22:25, 47:94]
    var_metric(:,k) = patn_neuropsych_tle.(patn_neuropsych_tle.Properties.VariableNames{i});
    var_name{k} = patn_neuropsych_tle.Properties.VariableNames{i};
    k = k+1;
end

SUBNO_anim_pt = [];
for i=1:length(sub_pt)
    SUBNO_anim_pt(i) = str2double(sub_pt{i}(3:end));
end

[TLE_sub,IA,IB] = intersect(SUBNO_anim_pt, SUBNO);
TLESide_sel = TLESide(IB);
AEDcount_sel = AEDcount(IB);

TLE_left = find(TLESide_sel == 'Left'); size(TLE_left)
TLE_right = find(TLESide_sel == 'Right');  size(TLE_right)
TLE_bilat = find(TLESide_sel == 'Bilateral');

disp(['tle, left:', num2str(length(TLE_left))])
disp(['tle, right:', num2str(length(TLE_right))])
% disp(['tle, bilat:', num2str(length(TLE_bilateral))])

TLE_idx = [];
TLE_idx.TLE_left = TLE_left;
TLE_idx.TLE_right = TLE_right;
TLE_idx.TLE_bilat = TLE_bilat;

% LI_anim_pt_val_tle =  LI_anim_pt_val(:,IA,:);
% LI_symb_pt_val_tle =  LI_symb_pt_val(:,IA,:);