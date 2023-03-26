
clc,

wi = []; w1 = 0; l = 0.1; ov = 0.01; j=1; %ov = l.*0.3
while w1+l < 2
    wi(j,:) = [w1, w1+l]; j=j+1; w1 = w1 + ov;
end
length(wi)

%%
data_save_dir = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/results/LI_subs';
cd(data_save_dir)
LI_anim_hc = load('LI_anim-hc');
LI_anim_pt = load('LI_anim-pt');

%%
% net_sel_mutiple_label = {'Angular'; 'Frontal'; 'Occipital'; 'Other'; 'PCingPrecun';'Temporal'};
net_sel_mutiple_label = {'Angular'; 'Frontal'; 'Occipital'; 'Other'; 'PCingPrecun';'Temporal'; 'BTLA'; 'VWFA'};

%%
LI_sub = LI_anim_hc.LI_sub; s_tag = 'anim-hc';
mLI_sub1 = squeeze(mean(LI_sub,2));

LI_sub = LI_anim_pt.LI_sub; s_tag = 'anim-pt';
mLI_sub2 = squeeze(mean(LI_sub,2));

%%
addpath('/data/MEG/Vahab/Github/MCW-MEGlab/FT/functions/External/brewermap');
network_sel = [1:3,6:8];
colr = distinguishable_colors(length(network_sel));

clc
% close all
mLI_sub = mLI_sub1 - mLI_sub2;

figure,
clear LI_val
for j=1:length(network_sel)
    hold on
    plot(mLI_sub(network_sel(j),:),'LineWidth',3, 'color', colr(j,:)),
    val = round(mean(wi(:,1),2),2);
    set(gca,'Xtick', 1:2:length(wi),'XtickLabel',val(1:2:end));
    set(gca,'FontSize',8,'XTickLabelRotation',90);
    set(gcf, 'Position', [1000   400   1100   500]);
end
% plot(mean(mLI_sub),'LineWidth',3),
lgnd = legend([net_sel_mutiple_label(network_sel); 'mean']);
title(['Diff LIs: HC - PT'])
ylabel('LI')
xlabel('time')
set(gca,'color','none');
set(lgnd,'color','none');

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