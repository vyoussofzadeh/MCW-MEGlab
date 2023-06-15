function sub_TLE_sub_data = ecpfunc_read_sub_TLE_sub(cfg)

sub_demog_data = cfg.sub_demog_data;
patn_neuropsych_data = cfg.patn_neuropsych_data;

%%
clear TLESide_sub_anim_pt
for i=1:length(sub_demog_data.sub_anim_pt)
    tmp = str2double(sub_demog_data.sub_anim_pt{i}(3:end));
    [C,IA,IB] = intersect(tmp,patn_neuropsych_data.SUBNO);
    TLESide_sub_anim_pt(i) = patn_neuropsych_data.TLESide(IB);
end

clear TLESide_sub_symb_patn
for i=1:length(sub_demog_data.sFiles_symb_patn)
    tmp = str2double(sub_demog_data.sub_symb_pt{i}(3:end));
    [C,IA,IB] = intersect(tmp,patn_neuropsych_data.SUBNO);
    TLESide_sub_symb_pt(i) = patn_neuropsych_data.TLESide(IB);
end

TLESide_sub_symb_patn_Left = find(contains(string(TLESide_sub_symb_pt), 'Left')==1);
TLESide_sub_symb_patn_Right = find(contains(string(TLESide_sub_symb_pt), 'Right')==1);

TLESide_sub_anim_patn_Left = find(contains(string(TLESide_sub_anim_pt), 'Left')==1);
TLESide_sub_anim_patn_Right = find(contains(string(TLESide_sub_anim_pt), 'Right')==1);

%%
sub_TLE_sub_data = [];
sub_TLE_sub_data.TLESide_sub_anim_pt = TLESide_sub_anim_pt;
sub_TLE_sub_data.TLESide_sub_symb_patn = TLESide_sub_symb_pt;

sub_TLE_sub_data.TLESide_sub_symb_patn_Left = TLESide_sub_symb_patn_Left;
sub_TLE_sub_data.TLESide_sub_symb_patn_Right = TLESide_sub_symb_patn_Right;
sub_TLE_sub_data.TLESide_sub_anim_patn_Left = TLESide_sub_anim_patn_Left;
sub_TLE_sub_data.TLESide_sub_anim_patn_Right = TLESide_sub_anim_patn_Right;


%%
% subjs_3 = cfg.subjs_3;
% subjs_2 = cfg.subjs_2;
% 
% sFiles_3 = cfg.sFiles_3;
% sFiles_2 = cfg.sFiles_2;
% 
% %%
% ECP_scriptdir = '/data/MEG/Research/ECP/Behavioural/processed';
% load(fullfile(ECP_scriptdir,'sub_demog.mat'));
% 
% k=1; ib_3 = [];
% for j=1:length(subjs_3)
%     [~, ~,ib] = intersect(subjs_3{j},sub_demog_save(:,1));
%     if ~isempty(ib)
%         ib_3(k) = ib;
%         k=k+1;
%     end
% end
% 
% k=1; ib_2 = [];
% for j=1:length(subjs_2)
%     [~, ~,ib] = intersect(subjs_2{j},sub_demog_save(:,1));
%     if ~isempty(ib)
%         ib_2(k) = ib;
%         k=k+1;
%     end
% end
% 
% sub_cond_3 = sub_cond_val(ib_3); sub_cond_2 = sub_cond_val(ib_2);
% 
% idx_ctrl_3 = find(sub_cond_3 ==1); idx_patn_3 = find(sub_cond_3 ==2);
% idx_ctrl_2 = find(sub_cond_2 ==1); idx_patn_2 = find(sub_cond_2 ==2);
% 
% sub_anim_hc = subjs_3(idx_ctrl_3);
% sub_symb_hc = subjs_2(idx_ctrl_2);
% sub_anim_pt = subjs_3(idx_patn_3);
% sub_symb_pt = subjs_2(idx_patn_2);
% 
% sFiles_anim_patn = sFiles_3(idx_patn_3); sFiles_symb_patn = sFiles_2(idx_patn_2);
% sFiles_anim_ctrl = sFiles_3(idx_ctrl_3); sFiles_symb_ctrl = sFiles_2(idx_ctrl_2);
% 
% 
% %%
% sub_demog_data = [];
% sub_demog_data.sFiles_anim_patn = sFiles_anim_patn;
% sub_demog_data.sFiles_anim_ctrl = sFiles_anim_ctrl;
% sub_demog_data.sFiles_symb_patn = sFiles_symb_patn;
% sub_demog_data.sFiles_symb_ctrl = sFiles_symb_ctrl;


% subjs_3 = cfg.subjs_3;
% subjs_2 = cfg.subjs_2;
% 
% sFiles_3 = cfg.sFiles_3;
% sFiles_2 = cfg.sFiles_2;
% 
% %%
% ECP_scriptdir = '/data/MEG/Research/ECP/Behavioural/processed';
% load(fullfile(ECP_scriptdir,'sub_demog.mat'));
% 
% k=1; ib_3 = [];
% for j=1:length(subjs_3)
%     [~, ~,ib] = intersect(subjs_3{j},sub_demog_save(:,1));
%     if ~isempty(ib)
%         ib_3(k) = ib;
%         k=k+1;
%     end
% end
% 
% k=1; ib_2 = [];
% for j=1:length(subjs_2)
%     [~, ~,ib] = intersect(subjs_2{j},sub_demog_save(:,1));
%     if ~isempty(ib)
%         ib_2(k) = ib;
%         k=k+1;
%     end
% end
% 
% sub_cond_3 = sub_cond_val(ib_3); sub_cond_2 = sub_cond_val(ib_2);
% 
% idx_ctrl_3 = find(sub_cond_3 ==1); idx_patn_3 = find(sub_cond_3 ==2);
% idx_ctrl_2 = find(sub_cond_2 ==1); idx_patn_2 = find(sub_cond_2 ==2);
% 
% sub_anim_hc = subjs_3(idx_ctrl_3);
% sub_symb_hc = subjs_2(idx_ctrl_2);
% sub_anim_pt = subjs_3(idx_patn_3);
% sub_symb_pt = subjs_2(idx_patn_2);
% 
% sFiles_anim_patn = sFiles_3(idx_patn_3); sFiles_symb_patn = sFiles_2(idx_patn_2);
% sFiles_anim_ctrl = sFiles_3(idx_ctrl_3); sFiles_symb_ctrl = sFiles_2(idx_ctrl_2);
% 
% 
% %%
% sub_demog_data = [];
% sub_demog_data.sFiles_anim_patn = sFiles_anim_patn;
% sub_demog_data.sFiles_anim_ctrl = sFiles_anim_ctrl;
% sub_demog_data.sFiles_symb_patn = sFiles_symb_patn;
% sub_demog_data.sFiles_symb_ctrl = sFiles_symb_ctrl;