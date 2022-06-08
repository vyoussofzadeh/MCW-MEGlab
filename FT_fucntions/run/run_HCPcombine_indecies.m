% for i=1:length(groups)
%     grois = groups{i};
%     idx = [];
%     for j=1:1
%         if contains(grois, '???')
%             idx(j) = strmatch(grois{j},lhctable.struct_names);
%         else
%             idx(j) = strmatch(['L_',grois{j}, '_ROI'],lhctable.struct_names);
%         end
%     end
%     idx_L{i} = idx;
% end
% idx_R = idx_L;

% bnd_L.label = table2cell(label_lr);
%

bnd_L.label = lhctable.struct_names;
bnd_R.label = rhctable.struct_names;

bnd_L.idx = lhix;
bnd_R.idx = rhix;

%%
% 
idx_L = []; idx_R = []; idx_clr_L = []; idx_clr_R = [];
k=1;
for i = 1:length(bnd_L.label)
    idx = [];
    if contains(bnd_L.label{i}, 'L_')
%         idx_label = find(strcmp(bnd_L.label, g_roi{i})==1); 
        idx = find(bnd_L.idx == i);
        idx_L = [idx_L;idx];
%         idx_clr_L = [idx_clr_L; repmat(colr_1(k,:),  length(idx), 1)];
        k=k+1;
    elseif contains(bnd_R.label{i}, 'R_')
%         idx_label = find(strcmp(bnd_R.label, g_roi{i})==1); 
        idx = find(bnd_R.idx == i);
        idx_R = [idx_R;idx];
%         idx_clr_R = [idx_clr_R; repmat(colr_1(k,:),  length(idx), 1)];
        k=k+1;
    end
end

%%

for i=1:length(groups)
    idx = find(lhix == i);
    mni_L(i,:) = mean(bnd_L.pnt(idx,:));
end

for i=1:length(groups)
    idx = find(rhix == i);
    mni_R(i,:) = mean(bnd_R.pnt(idx,:));
end

mni = [];
for i=1:length(groups)
    mni = [mni;mni_L(i,:)];
    mni = [mni;mni_R(i,:)];
end

if flag.plotatlas ==1
    % close all
    figure;
    ft_plot_mesh(bnd_L, 'faecolor', 'brain',  'vertexcolor', ...
        lhix, 'facealpha', 1);
    hold on,
    ft_plot_mesh(bnd_R, 'faecolor', 'brain',  'vertexcolor', ...
        rhix, 'facealpha', 1);
    title('HCP-MMP1 - combined(362 ROIs)')
end