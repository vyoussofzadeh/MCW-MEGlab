function sub_TLE_sub_data = ecpfunc_read_sub_TLE_sub_contrast(cfg)

sub_demog_data = cfg.sub_demog_data;
patn_neuropsych_data = cfg.patn_neuropsych_data;

%%
clear TLESide_sub_anim_pt
for i=1:length(sub_demog_data.sub_pt)
    tmp = str2double(sub_demog_data.sub_pt{i}(3:end));
    [C,IA,IB] = intersect(tmp,patn_neuropsych_data.SUBNO);
    TLESide_sub_pt(i) = patn_neuropsych_data.TLESide(IB);
end

TLESide_sub_patn_Left = find(contains(string(TLESide_sub_pt), 'Left')==1);
TLESide_sub_patn_Right = find(contains(string(TLESide_sub_pt), 'Right')==1);

%%
sub_TLE_sub_data = [];
sub_TLE_sub_data.TLESide_sub_pt = TLESide_sub_pt;

sub_TLE_sub_data.TLESide_sub_patn_Left = TLESide_sub_patn_Left;
sub_TLE_sub_data.TLESide_sub_patn_Right = TLESide_sub_patn_Right;
