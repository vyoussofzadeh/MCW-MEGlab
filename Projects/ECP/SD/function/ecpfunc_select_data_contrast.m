function S_data = ecpfunc_select_data_contrast(cfg)


% taskcon = {'Anim', 'Symbol'};


switch cfg.select_data
    case 1
        J = 1; subcon = 'Ctrl'; sFiles_in = cfg.sub_demog_data.sFiles_ctrl; 
        sFiles_subid = cfg.sub_demog_data.sub_hc;  %- ctrl, Anim
    case 2
        J = 1; subcon = 'Patn'; sFiles_in = cfg.sub_demog_data.sFiles_patn; 
        sFiles_subid = cfg.sub_demog_data.sub_pt; %- Patn, Anim
    case 3
        J = 2; subcon = 'Ctrl'; sFiles_in = cfg.sub_demog_data.sFiles_ctrl; 
        sFiles_subid = cfg.sub_demog_data.sub_hc; %- Ctrl, Symb
    case 4
        J = 2; subcon = 'Patn'; sFiles_in = cfg.sub_demog_data.sFiles_patn; 
        sFiles_subid = cfg.sub_demog_data.sub_pt; %- Patn, Symb
end

% taskcon = {'Anim', 'Symbol'};
% 
% sFiles_anim_hc = 'Group_analysis/1_LCMV_Subjects/results_average_230110_2174.mat';
% sFiles_anim_pt = 'Group_analysis/1_LCMV_Subjects/results_average_230110_2175.mat';
% sFiles_symb_hc = 'Group_analysis/1_LCMV_Subjects/results_average_230111_1409.mat';
% sFiles_symb_pt = 'Group_analysis/1_LCMV_Subjects/results_average_230113_1150.mat';
% % sFiles_anim_hc_sFiles_anim_pt = 'Group_analysis/@intra/results_221226_1142.mat';
% 
% % disp('1: Anim, Ctrl')
% % disp('2: Anim, Patn')
% % disp('3: Symbol, Ctrl')
% % disp('4: Symbol, Patn')
% switch cfg.select_data
%     case 1
%         s_avg_Input = sFiles_anim_hc; s_tag = 'anim-hc';
%     case 2
%         s_avg_Input = sFiles_anim_pt; s_tag = 'anim-pt';
%     case 3
%         s_avg_Input = sFiles_symb_hc; s_tag = 'symb-hc';
%     case 4
%         s_avg_Input = sFiles_symb_pt; s_tag = 'symb-pt';
% end

S_data = [];
S_data.sFiles_in = sFiles_in;
% S_data.s_avg_Input = s_avg_Input;
% S_data.s_tag = s_tag;
S_data.sFiles_subid = sFiles_subid;
S_data.subcon = subcon;
% S_data.taskcon = taskcon{J};


end