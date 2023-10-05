
subdir = '/data/MEG/Vahab/Github/MCW-MEGlab/Projects/ECP/behavioural_demog/Subj_demog';
cd(subdir)
%% reading subject demographic information
clc
[num, txt, sub_demog] = xlsread(fullfile(subdir,'ECP Progress Tracker'));

sub_list = []; k=1;
sub_cond = [];
sub_cond_idx = [];
for i=1:length(sub_demog)
    tmp = sub_demog(i,3);
    tmp2 = sub_demog(i,9);
    if ~isnan(tmp{1})
        if contains(tmp{1}, 'EC') && ~sum(isnan(tmp2{1}))
            sub_list{k} = tmp{1};
            sub_cond{k} = tmp2{1};
            if contains(sub_cond{k},'Ctrl')
                sub_cond_idx{k} = 1; sub_cond_val(k) = 1;
            elseif contains(sub_cond{k},'Patn')
                sub_cond_idx{k} = 2; sub_cond_val(k) = 2;
            end
            
            k=k+1;
        end
    end
end
disp(sub_list') % Suject ID
disp(sub_cond') % Sex
% disp(sub_cond_idx')

idx_ctrl = find(sub_cond_val ==1);
idx_patn = find(sub_cond_val ==2);

sub_demog_save = [sub_list',sub_cond',sub_cond_idx'];
save(fullfile(subdir,'sub_demog.mat'),'sub_demog_save','sub_cond_val')

%% Handedness
clear sub_list sub_sex sub_demog_ctrl_save sub_EHQ 

clc
[num, txt, sub_demog] = xlsread(fullfile(subdir,'Demographics_control_20200626'));

sub_list = []; k=1;
sub_cond = [];
sub_cond_idx = [];
for i=1:length(sub_demog)
    tmp = sub_demog(i,1); tmp2 = sub_demog(i,2); tmp3 = sub_demog(i,3); tmp4 = sub_demog(i,15);
    if ~isnan(tmp{1})
        if contains(tmp{1}, 'EC') && ~sum(isnan(tmp2{1}))
            sub_list{k} = tmp{1};
            sub_sex{k} = tmp2{1};
            sub_age{k} = tmp3{1};
            sub_EHQ{k} = tmp4{1};
            k=k+1;
        end
    end
end
disp(sub_sex')
disp(sub_list')
disp(sub_age')
% disp(sub_cond_idx')

sub_demog_ctrl_save = [sub_list',sub_sex', sub_age', sub_EHQ'];

save(fullfile(subdir,'sub_demog_ctrl.mat'),'sub_demog_ctrl_save')

%%
clear sub_list sub_sex sub_TLEside sub_demog_PT_save sub_EHQ

clc
[num, txt, sub_demog] = xlsread(fullfile(subdir,'Demographics_TLE_20200626'));

sub_list = []; k=1;
sub_cond = [];
sub_cond_idx = [];
for i=1:length(sub_demog)
    tmp = sub_demog(i,1); tmp1 = sub_demog(i,2); tmp2 = sub_demog(i,3); tmp3 = sub_demog(i,4); tmp4 = sub_demog(i,16);
    if ~isnan(tmp{1})
        if contains(tmp{1}, 'EC') && ~sum(isnan(tmp1{1}))
            sub_list{k} = tmp{1};
            sub_TLEside{k} = tmp1{1};
            sub_sex{k} = tmp2{1};
            sub_age{k} = tmp3{1};
            sub_EHQ{k} = tmp4{1};
            k=k+1;
        end
    end
end
disp(sub_sex')
disp(sub_list')
disp(sub_TLEside')
disp(sub_age')
% disp(sub_cond_idx')

sub_demog_PT_save = [sub_list',sub_sex',sub_TLEside',sub_age,sub_EHQ'];

save(fullfile(subdir,'sub_demog_pnts.mat'),'sub_demog_PT_save');


%%
clear sub_list sub_sex sub_TLEside sub_demog_PT_save sub_EHQ sub_WASI_VocR

clc
[num, txt, sub_demog] = xlsread(fullfile(subdir,'ECP_FINAL Data Set_v_04_04_19'));

sub_list = []; k=1;
sub_cond = [];
sub_cond_idx = [];
for i=1:length(sub_demog)
    tmp = sub_demog(i,1); tmp1 = sub_demog(i,119); tmp2 = sub_demog(i,206); tmp3 = sub_demog(i,207); tmp4 = sub_demog(i,208); tmp5 = sub_demog(i,68); 
    if ~isnan(tmp{1})
        if contains(tmp{1}, 'EC') && ~sum(isnan(tmp1{1}))
            sub_list{k}          = tmp{1};
            sub_BNT_T{k}         = tmp1{1};
            sub_VOCABagecorSS{k} = tmp2{1};
            sub_VOCABtscore{k}   = tmp3{1};
            sub_ORALRagecorSS{k} = tmp4{1};
            sub_WASI_VocR{k}     = tmp5{1};
            k=k+1;
        end
    end
end
disp(sub_list')
disp(sub_BNT_T')
disp(sub_VOCABagecorSS')
disp(sub_VOCABtscore')
disp(sub_WASI_VocR')
% disp(sub_cond_idx')

sub_behave_save = [sub_list',sub_BNT_T',sub_VOCABagecorSS',sub_VOCABtscore',sub_ORALRagecorSS',sub_WASI_VocR'];

save(fullfile(subdir,'sub_behave.mat'),'sub_behave_save');

