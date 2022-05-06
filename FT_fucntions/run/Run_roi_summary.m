
%%
Run_aal_labels

clear label
for i=1:length(idx2)
    label{i}  = aal_label1{idx2(i)};
end

clear DataArray2
DataArray2 = DataArray1(idx2,:);

figure,
bar_handle = bar(DataArray2,'grouped','BarWidth', 1);
% set(bar_handle, 'edgecolor','none')
set(bar_handle(1),'FaceColor',Colors(1,:))
set(bar_handle(2),'FaceColor',Colors(2,:))
L = length(DataArray2);
% set(gca,'Xtick', 1:L,'XtickLabel',label); set(gca,'FontSize',10,'XTickLabelRotation',45)
set(gca,'Xtick', 1:L,'XtickLabel',[1:L/2,1:L/2]); set(gca,'FontSize',10,'XTickLabelRotation',90)
grid off
box off
set(gca,'color','none');
% axis on
xlim([0,L+1])
set(gca,'FontName','HelveticaNeueLT Std Lt');
xlabel(xtag,'FontSize', 12);
ylabel(ytag,'FontSize', 12);
set(gcf, 'Position', [1000   500   1000  400]);
% title('Group average')

mu = max(DataArray1(:));
hline = refline([0 0.8.*mu]);
hline.Color = 'r';
lgnd = legend({'Aud Def Naming','Vis Pic Naming','Threshold'});
set(lgnd,'color','none');


%% sort
% [val,idx] = sort((m),'descend');

% if abs(min(m)) > abs(max(m))
%     m1 = -m;
%     idx2 = find(m1 > thre.*max(m1));
% else
%     idx2 = find(m > thre.*max(m));
% end
% n = length(idx2);
% col = ones(n,3);
% col(:,[1,3]) = 0;
% for i=1:n
%     m_sig = zeros(size(m,1),size(m,2));
%     m_sig(idx2(i)) = m(idx2(i));
%     c = bar(m_sig);
%     set(c,'faceColor',col(i,:));
%     set(c,'BarWidth',1);
% end
% title(['ROIs (',num2str(100*thre),'% threshold)']);
% 
% %% roi summary
% L = length(m);
% id = table([1:L]','VariableNames',{'ID'});
% Z = table(m,'VariableNames',{mask});
% ROI = [id,cell2table(label),Z];
% coor_var = table(coor,'VariableNames',{'MNI'});
% ROI(:,4) = coor_var;
% ROI.Properties.VariableNames{'Var4'} = 'MNI';
%

%% roi selected
% [~,idx] = sort((abs(m(idx2))),'descend');
% % disp(ROI)
% ROI_sel = ROI(idx2(idx),:);
% disp(ROI_sel)