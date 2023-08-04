clc, clear

wi = []; w1 = 0; l = 0.1; ov = 0.01; j=1; %ov = l.*0.3
while w1+l < 2
    wi(j,:) = [w1, w1+l]; j=j+1; w1 = w1 + ov;
end
length(wi)

%%
% net_sel_mutiple_label = {'Angular'; 'Frontal'; 'Occipital'; 'Other'; 'PCingPrecun';'Temporal'};
net_sel_mutiple_label = {'Angular'; 'Frontal'; 'Occipital'; 'Other'; 'PCingPrecun';'Temporal'; 'BTLA'; 'VWFA'};
% net_sel_mutiple_label = {'Angular'; 'Frontal'; 'Temporal'};


%%
network_sel = [1:3,6:8];
colr = distinguishable_colors(length(network_sel));

%%
data_save_dir = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/results/LI_subs/group_8net_100ms';
cd(data_save_dir)

%%
cd(data_save_dir)

LI_anim_hc = load('LI_anim-hc');
LI_anim_pt = load('LI_anim-pt');

LI_symb_hc = load('LI_symb-hc');
LI_symb_pt = load('LI_symb-pt');

%% anim vs. symb
close all

mLI_sub1 = squeeze(mean(LI_anim_hc.LI_sub,2));
mLI_sub2 = squeeze(mean(LI_symb_hc.LI_sub,2)); tag = 'anim vs. symb, hc';

clc
mLI_sub_hc = mLI_sub1 - mLI_sub2;

figure,
subplot(3,1,1)
clear LI_val
for j=1:length(network_sel)
    hold on
    plot(mLI_sub_hc(network_sel(j),:),'LineWidth',3, 'color', colr(j,:)),
    val = round(mean(wi(:,1),2),2);
    set(gca,'Xtick', 1:2:length(wi),'XtickLabel',val(1:2:end));
    set(gca,'FontSize',8,'XTickLabelRotation',90);
    set(gcf, 'Position', [1000   400   1100   500]);
end
lgnd = legend([net_sel_mutiple_label(network_sel); 'mean']);
title(tag)
ylabel('LI')
xlabel('time')
set(gca,'color','none');
set(lgnd,'color','none');

mLI_sub1 = squeeze(mean(LI_anim_pt.LI_sub,2));
mLI_sub2 = squeeze(mean(LI_symb_pt.LI_sub,2)); tag = 'anim vs. symb, pt';

mLI_sub_pt = mLI_sub1 - mLI_sub2;

% figure,
subplot(3,1,2)
clear LI_val
for j=1:length(network_sel)
    hold on
    plot(mLI_sub_pt(network_sel(j),:),'LineWidth',3, 'color', colr(j,:)),
    val = round(mean(wi(:,1),2),2);
    set(gca,'Xtick', 1:2:length(wi),'XtickLabel',val(1:2:end));
    set(gca,'FontSize',8,'XTickLabelRotation',90);
end
lgnd = legend([net_sel_mutiple_label(network_sel); 'mean']);
title(tag)
ylabel('LI')
xlabel('time')
set(gca,'color','none');
set(lgnd,'color','none');


mLI_sub_diff = mLI_sub_hc - mLI_sub_pt; tag = 'hc - pt';
%% LI ftp (fronto-tempro-pari network)
% % % close all
% % 
% % LI_ftp_anim_hc = load('LI-ftp-anim-hc');
% % LI_ftp_anim_pt = load('LI-ftp-anim-pt');
% % 
% % LI_ftp_symb_hc = load('LI-ftp-symb-hc');
% % LI_ftp_symb_pt = load('LI-ftp-symb-pt');
% % 
% % cfg = [];
% % cfg.sub_sel = LI_anim_hc.sFiles_subid;
% % cfg.d_in = LI_ftp_anim_hc.m_LI_sub - LI_ftp_symb_hc.m_LI_sub;
% % cfg.tit = 'anim vs. symb, hc';
% % do_barplot_LI(cfg)
% % set(gcf, 'Position', [1000   400   500   300]);
% % disp(cfg.d_in)
% % 
% % [C,IA,IB] = intersect(LI_anim_pt.sFiles_subid, LI_symb_pt.sFiles_subid);
% % 
% % cfg = [];
% % cfg.sub_sel = LI_anim_pt.sFiles_subid(IA);
% % cfg.d_in = LI_ftp_anim_pt.m_LI_sub(IA)- LI_ftp_symb_pt.m_LI_sub(IB);
% % % cfg.d_in = LI_ftp_anim_pt.m_LI_sub(IA);
% % cfg.tit = 'anim vs. symb, pt';
% % do_barplot_LI(cfg)
% % set(gcf, 'Position', [1000   400   1000   300]);

%%
% clc, close all
% subid = 20;
% 
% cfg = [];
% cfg.net_sel_mutiple_label = net_sel_mutiple_label;
% cfg.LI_sub = LI_anim_hc.LI_sub - LI_symb_hc.LI_sub;
% cfg.wi = wi;
% cfg.subsel = subid;
% mm = do_plot_sub_LI(cfg);
% disp(mean(mm(1:100)))


%%
%
% % net_sel_mutiple_label = LI_anim_hc.net_sel_mutiple_label;
%
% mLI_sub = mLI_sub1 - mLI_sub2;
%
% colr = distinguishable_colors(length(net_sel_mutiple_label));
%
% clc
% % close all
% % mLI_sub = squeeze(mean(LI_sub,2));
% figure,
% for j=1:length(net_sel_mutiple_label)
%     hold on
%     plot(mLI_sub(j,:),'LineWidth',3, 'color', colr(j,:)),
%     val = round(mean(wi(:,1),2),2);
%     set(gca,'Xtick', 1:2:length(wi),'XtickLabel',val(1:2:end));
%     set(gca,'FontSize',8,'XTickLabelRotation',90);
%     set(gcf, 'Position', [1000   400   1100   500]);
% %     set(gca,'color','none');
%
% end
% % plot(mean(mLI_sub),'LineWidth',3),
% lgnd = legend([net_sel_mutiple_label; 'mean'])
% title(['mean subject LIs: ', s_tag])
% % lgnd = legend(net_sel_mutiple_label(network_sel));
% % title(['subject LIs: ', s_tag,])
% set(gca,'color','none');
% set(lgnd,'color','none');
%
% %%
% figure,
% for j=1:size(mm,1)
%     hold on
%     plot(mm(j,:),'LineWidth',3, 'color', colr(j,:)),
%     val = round(mean(wi(:,1),2),2);
%     set(gca,'Xtick', 1:2:length(wi),'XtickLabel',val(1:2:end));
%     set(gca,'FontSize',8,'XTickLabelRotation',90);
%     set(gcf, 'Position', [1000   400   1100   500]);
% end
% lgnd = legend(net_sel_mutiple_label(network_sel));
% title(['subject LIs: ', s_tag,])
% set(gca,'color','none');
% set(lgnd,'color','none');
% figure,
subplot(3,1,3)
clear LI_val
for j=1:length(network_sel)
    hold on
    plot(mLI_sub_diff(network_sel(j),:),'LineWidth',3, 'color', colr(j,:)),
    val = round(mean(wi(:,1),2),2);
    set(gca,'Xtick', 1:2:length(wi),'XtickLabel',val(1:2:end));
    set(gca,'FontSize',8,'XTickLabelRotation',90);
end
lgnd = legend([net_sel_mutiple_label(network_sel); 'mean']);
title(tag)
ylabel('LI')
xlabel('time')
set(gca,'color','none');
set(lgnd,'color','none');

set(gcf, 'Position', [1000   400   1000   900]);

%%
mLI_sub_mdiff = mean(mLI_sub_hc,1) - mean(mLI_sub_pt,1); tag = 'hc - pt, mean';

figure,
plot(mLI_sub_mdiff,'LineWidth',3, 'color','k'),
val = round(mean(wi(:,1),2),2);
set(gca,'Xtick', 1:2:length(wi),'XtickLabel',val(1:2:end));
set(gca,'FontSize',8,'XTickLabelRotation',90);
set(gcf, 'Position', [1000   400   1100   300]);
title(tag)
set(gca,'color','none');
set(lgnd,'color','none');

%%
patn_neuropsych_data = ecpfunc_read_patn_neuropsych();


