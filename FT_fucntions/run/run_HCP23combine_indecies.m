for i=1:length(groups)
    grois = groups{i}(2:end);
    idx = [];
    for j=1:length(grois)
        if contains(grois, '???')
            idx(j) = strmatch(grois{j},lhctable.struct_names);
        else
            idx(j) = strmatch(['L_',grois{j}, '_ROI'],lhctable.struct_names);
        end
    end
    idx_L{i} = idx;
end
idx_R = idx_L;

lhvtx_combined = lhvtx;
rhvtx_combined = rhvtx;

lhix_combined = lhix;
rhix_combined = rhix;

% colr_combined = hsv(length(idx_L));
lhix_combined = zeros(size(lhix,1),1);
rhix_combined = zeros(size(rhix,1),1);

for i=1:length(idx_L)
    for j=1:length(idx_L{i})
        idx = find(lhix == idx_L{i}(j));
        lhix_combined(idx,:) = i;
    end
end

for i=1:length(idx_R)
    for j=1:length(idx_R{i})
        idx = find(rhix == idx_R{i}(j));
        rhix_combined(idx,:) = i;
    end
end

for i=1:length(idx_L)
    idx = find(lhix_combined == i);
    mni_L(i,:) = mean(bnd_L.pnt(idx,:));
end

for i=1:length(idx_R)
    idx = find(rhix_combined == i);
    mni_R(i,:) = mean(bnd_R.pnt(idx,:));
end

mni = [];
for i=1:length(groups)
    mni = [mni;mni_L(i,:)];
    mni = [mni;mni_R(i,:)];
end

bnd_L.idx_23 = lhix_combined;
bnd_R.idx_23 = rhix_combined;


if flag.plotatlas ==1
    % close all
    figure;
    ft_plot_mesh(bnd_L, 'faecolor', 'brain',  'vertexcolor', ...
        lhix_combined, 'facealpha', 1);
    hold on,
    ft_plot_mesh(bnd_R, 'faecolor', 'brain',  'vertexcolor', ...
        rhix_combined, 'facealpha', 1);
    title('HCP-MMP1 - combined(23 ROIs)')
end