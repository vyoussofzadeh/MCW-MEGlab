function [LI_tle, Pow_tle] = do_update_tle_li_pow(cfg_main)


LI = cfg_main.LI;
tle_IA = cfg_main.tle_IA;
TLE_idx = cfg_main.TLE_idx;

%%
LI_anim_pt_val_tle =  LI.LI_anim_pt.LI_sub(:,tle_IA,:);
LI_symb_pt_val_tle =  LI.LI_symb_pt.LI_sub(:,tle_IA,:);

% tle left
LI_anim_pt_val_tle_left = LI_anim_pt_val_tle(:,TLE_idx.TLE_left,:);
LI_symb_pt_val_tle_left = LI_symb_pt_val_tle(:,TLE_idx.TLE_left,:);

% tle right
LI_anim_pt_val_tle_right = LI_anim_pt_val_tle(:,TLE_idx.TLE_right,:);
LI_symb_pt_val_tle_right = LI_symb_pt_val_tle(:,TLE_idx.TLE_right,:);

% tle bilateral
LI_anim_pt_val_tle_bilat = LI_anim_pt_val_tle(:,TLE_idx.TLE_bilat,:);
LI_symb_pt_val_tle_bilat = LI_symb_pt_val_tle(:,TLE_idx.TLE_bilat,:);

%% Power
pow_anim_pt_val_LH_tle = LI.LI_anim_pt.pow_LH(:,tle_IA,:);
pow_symb_pt_val_LH_tle = LI.LI_symb_pt.pow_LH(:,tle_IA,:);

% tle left
pow_anim_pt_val_LH_tle_left = pow_anim_pt_val_LH_tle(:,TLE_idx.TLE_left,:);
pow_symb_pt_val_LH_tle_left = pow_symb_pt_val_LH_tle(:,TLE_idx.TLE_left,:);

% tle right
pow_anim_pt_val_LH_tle_right = pow_anim_pt_val_LH_tle(:,TLE_idx.TLE_right,:);
pow_symb_pt_val_LH_tle_right = pow_symb_pt_val_LH_tle(:,TLE_idx.TLE_right,:);

% tle bilateral
pow_anim_pt_val_LH_tle_bilat = pow_anim_pt_val_LH_tle(:,TLE_idx.TLE_bilat,:);
pow_symb_pt_val_LH_tle_bilat = pow_symb_pt_val_LH_tle(:,TLE_idx.TLE_bilat,:);

%%
LI_tle = LI;

LI_tle.LI_anim_pt.LI_sub_tle = LI_anim_pt_val_tle;
LI_tle.LI_symb_pt.LI_sub_tle = LI_symb_pt_val_tle;

LI_tle.LI_anim_pt.LI_sub_tle_left = LI_anim_pt_val_tle_left;
LI_tle.LI_symb_pt.LI_sub_tle_left = LI_symb_pt_val_tle_left;

LI_tle.LI_anim_pt.LI_sub_tle_right = LI_anim_pt_val_tle_right;
LI_tle.LI_symb_pt.LI_sub_tle_right = LI_symb_pt_val_tle_right;

LI_tle.LI_anim_pt.LI_sub_tle_bilat = LI_anim_pt_val_tle_bilat;
LI_tle.LI_symb_pt.LI_sub_tle_bilat = LI_symb_pt_val_tle_bilat;

%%
LI_tle.LI_anim_pt.pow_sub_tle_LH = pow_anim_pt_val_LH_tle;
LI_tle.LI_symb_pt.pow_sub_tle_LH = pow_symb_pt_val_LH_tle;

LI_tle.LI_anim_pt.pow_sub_tle_LH_left = pow_anim_pt_val_LH_tle_left;
LI_tle.LI_symb_pt.pow_sub_tle_LH_left = pow_symb_pt_val_LH_tle_left;

LI_tle.LI_anim_pt.pow_sub_tle_LH_right = pow_anim_pt_val_LH_tle_right;
LI_tle.LI_symb_pt.pow_sub_tle_LH_right = pow_symb_pt_val_LH_tle_right;

LI_tle.LI_anim_pt.pow_sub_tle_LH_bilat = pow_anim_pt_val_LH_tle_bilat;
LI_tle.LI_symb_pt.pow_sub_tle_LH_bilat = pow_symb_pt_val_LH_tle_bilat;

%%
% % Pow_tle = [];
% Pow_tle.pow_anim_pt_val_LH_tle_left = pow_anim_pt_val_LH_tle_left;
% Pow_tle.pow_symb_pt_val_LH_tle_left = pow_symb_pt_val_LH_tle_left;
% 
% Pow_tle.pow_anim_pt_val_LH_tle_right = pow_anim_pt_val_LH_tle_right;
% Pow_tle.pow_symb_pt_val_LH_tle_right = pow_symb_pt_val_LH_tle_right;
% 
% Pow_tle.pow_anim_pt_val_LH_tle_bilat = pow_anim_pt_val_LH_tle_bilat;
% Pow_tle.pow_symb_pt_val_LH_tle_bilat = pow_symb_pt_val_LH_tle_bilat;



