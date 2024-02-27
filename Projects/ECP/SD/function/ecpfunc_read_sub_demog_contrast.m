function sub_demog_data = ecpfunc_read_sub_demog_contrast(cfg)


subjs = cfg.subjs;

sFiles = cfg.sFiles;

%%
ECP_scriptdir = '/data/MEG/Research/ECP/Behavioural/processed';
load(fullfile(ECP_scriptdir,'sub_demog.mat'));

k=1; ib = [];
for j=1:length(subjs)
    [~, ~,ib] = intersect(subjs{j},sub_demog_save(:,1));
    if ~isempty(ib)
        ib_new(k) = ib;
%         disp(ib);
%         pause,
        k=k+1;
    end
end

%%
sub_cond = sub_cond_val(ib_new); 

idx_ctrl = find(sub_cond ==1); idx_patn = find(sub_cond ==2);

sub_hc = subjs(idx_ctrl);
sub_pt = subjs(idx_patn);

sFiles_patn = sFiles(idx_patn); 
sFiles_ctrl = sFiles(idx_ctrl); 


%%
sub_demog_data = [];
sub_demog_data.sFiles_patn = sFiles_patn;
sub_demog_data.sFiles_ctrl = sFiles_ctrl;
% sub_demog_data.sFiles_symb_patn = sFiles_symb_patn;
% sub_demog_data.sFiles_symb_ctrl = sFiles_symb_ctrl;
sub_demog_data.sub_hc = sub_hc;
% sub_demog_data.sub_symb_hc = sub_symb_hc;
sub_demog_data.sub_pt = sub_pt;
% sub_demog_data.sub_symb_pt = sub_symb_pt;