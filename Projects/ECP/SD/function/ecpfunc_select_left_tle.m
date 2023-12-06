function LI_pt_left_tle = ecpfunc_select_left_tle(cfg_main)


%% all pts
patn_neuropsych_tle = cfg_main.patn_neuropsych_tle;
LI_pt = cfg_main.LI_pt;
pt_ID = LI_pt.pt_ID;
LI_methods = LI_pt.LI_methods;

TLESide = patn_neuropsych_tle.TLESide; SUBNO = patn_neuropsych_tle.SUBNO;

SUBNO_pt = [];
for i=1:length(pt_ID)
    SUBNO_pt(i) = str2double(pt_ID{i}(3:end));
end

[~,IA,IB] = intersect(SUBNO_pt, SUBNO); % PT with TLE known side
TLESide_sel = TLESide(IB); % TLE side of known PT

%% tle pts
LI_pt_tle = LI_pt;

for i=1:length(LI_methods)
    LI_pt_val = LI_pt.(LI_methods{i}).LI_sub;
    LI_pt_val_max = LI_pt.(LI_methods{i}).m_LI_max_sub;
    
    LI_pt_tle.(LI_methods{i}).LI_sub = LI_pt_val(:,IA,:);
    LI_pt_tle.(LI_methods{i}).m_LI_max_sub = LI_pt_val_max(:,IA,:);
end

LI_pt_tle.pt_ID = LI_pt_tle.pt_ID(IA);

%% left tle pt
TLE_left = TLESide_sel == 'Left'; % Left_TLE idx

pt_ID_ltle = LI_pt_tle.pt_ID(TLE_left);


[~,~,IB] = intersect(pt_ID_ltle, LI_pt_tle.pt_ID);

%% Update left tle Lis 
LI_pt_left_tle = LI_pt_tle;

for i=1:length(LI_methods)
    LI_pt_val = LI_pt_left_tle.(LI_methods{i}).LI_sub;
    LI_pt_val_max = LI_pt_left_tle.(LI_methods{i}).m_LI_max_sub;
    
    LI_pt_left_tle.(LI_methods{i}).LI_sub = LI_pt_val(:,IB,:);
    LI_pt_left_tle.(LI_methods{i}).m_LI_max_sub = LI_pt_val_max(:,IB,:);
end

LI_pt_left_tle.pt_ID = LI_pt_tle.pt_ID(IB);

end




