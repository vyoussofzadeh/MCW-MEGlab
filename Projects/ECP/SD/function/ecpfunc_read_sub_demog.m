function sub_demog_data = ecpfunc_read_sub_demog(cfg)


subjs_3 = cfg.subjs_3;
subjs_2 = cfg.subjs_2;

sFiles_3 = cfg.sFiles_3;
sFiles_2 = cfg.sFiles_2;

%%
ECP_scriptdir = '/data/MEG/Research/ECP/Behavioural/processed';
load(fullfile(ECP_scriptdir,'sub_demog.mat'));

k=1; ib_3 = [];
for j=1:length(subjs_3)
    [~, ~,ib] = intersect(subjs_3{j},sub_demog_save(:,1));
    if ~isempty(ib)
        ib_3(k) = ib;
        k=k+1;
    end
end

k=1; ib_2 = [];
for j=1:length(subjs_2)
    [~, ~,ib] = intersect(subjs_2{j},sub_demog_save(:,1));
    if ~isempty(ib)
        ib_2(k) = ib;
        k=k+1;
    end
end

sub_cond_3 = sub_cond_val(ib_3); sub_cond_2 = sub_cond_val(ib_2);

idx_ctrl_3 = find(sub_cond_3 ==1); idx_patn_3 = find(sub_cond_3 ==2);
idx_ctrl_2 = find(sub_cond_2 ==1); idx_patn_2 = find(sub_cond_2 ==2);

sub_anim_hc = subjs_3(idx_ctrl_3);
sub_symb_hc = subjs_2(idx_ctrl_2);
sub_anim_pt = subjs_3(idx_patn_3);
sub_symb_pt = subjs_2(idx_patn_2);

sFiles_anim_patn = sFiles_3(idx_patn_3); sFiles_symb_patn = sFiles_2(idx_patn_2);
sFiles_anim_ctrl = sFiles_3(idx_ctrl_3); sFiles_symb_ctrl = sFiles_2(idx_ctrl_2);


%%
sub_demog_data = [];
sub_demog_data.sFiles_anim_patn = sFiles_anim_patn;
sub_demog_data.sFiles_anim_ctrl = sFiles_anim_ctrl;
sub_demog_data.sFiles_symb_patn = sFiles_symb_patn;
sub_demog_data.sFiles_symb_ctrl = sFiles_symb_ctrl;
sub_demog_data.sub_anim_hc = sub_anim_hc;
sub_demog_data.sub_symb_hc = sub_symb_hc;
sub_demog_data.sub_anim_pt = sub_anim_pt;
sub_demog_data.sub_symb_pt = sub_symb_pt;