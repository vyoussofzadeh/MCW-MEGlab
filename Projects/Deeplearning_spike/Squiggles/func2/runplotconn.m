% --- 1. Get layout for Neuromag306 sensors in YOUR channel order ---
cfg = [];
cfg.grad    = MEG_data.grad;      % neuromag306
cfg.channel = MEG_data.label;     % keep same order as outsum
lay = ft_prepare_layout(cfg);

pos = lay.pos;          % [306 x 2] 2D positions (x = left/right, y = ant/post)
labels = lay.label;     % should match MEG_data.label (plus maybe COMNT/SCALE)
nChan = numel(MEG_data.label);

% keep only real MEG channels (skip COMNT/SCALE if present)
labels = labels(1:nChan);
pos    = pos(1:nChan,:);

x = pos(:,1);
y = pos(:,2);

% --- 2. Define left/right and anterior/posterior stripes ---
xmid = median(x);
yq   = quantile(y,[0.25 0.5 0.75]);   % 4 stripes back?front

% ROI order: 1=L-Occ, 2=R-Occ, 3=L-Temp, 4=R-Temp, 5=L-CentPar,
%            6=R-CentPar, 7=L-Front, 8=R-Front
roi_names = {...
    'L-Occ', 'R-Occ', ...
    'L-Temp','R-Temp', ...
    'L-CentPar','R-CentPar', ...
    'L-Front','R-Front'};

roi_idx = zeros(nChan,1);  % ROI index per channel

for ch = 1:nChan
    % stripe in y (posterior?anterior)
    if y(ch) < yq(1)
        stripe = 1;   % occipital
    elseif y(ch) < yq(2)
        stripe = 2;   % temporal
    elseif y(ch) < yq(3)
        stripe = 3;   % centralparietal
    else
        stripe = 4;   % frontal
    end

    % hemisphere in x
    if x(ch) < xmid
        hemi = 1;     % left
    else
        hemi = 2;     % right
    end

    roi_idx(ch) = (stripe-1)*2 + hemi;  % 1..8
end

%%
% outsum: [306 x 306], same channel order as MEG_data.label

[roi_idx_sorted, sort_idx] = sort(roi_idx);
out_ord = outsum(sort_idx, sort_idx);

figure;
imagesc(out_ord);
axis square;
colorbar;
title('Connectivity grouped by Neuromag306 ROIs');

% group boundaries + centers
edges   = [find(diff(roi_idx_sorted)~=0); numel(roi_idx_sorted)];
starts  = [1; edges(1:end-1)+1];
centers = floor((starts + edges)/2);

xticks(centers);
yticks(centers);
xticklabels(roi_names);
yticklabels(roi_names);
xtickangle(45);
set(gca,'TickDir','out','FontSize',10);

hold on;
for e = edges(1:end-1)'
    xline(e+0.5,'k-');
    yline(e+0.5,'k-');
end
hold off;

%%

% outsum: [306 x nTime]
[roi_idx_sorted, sort_idx] = sort(roi_idx);
out_ord = outsum(sort_idx, :);

figure;
imagesc(MEG_data.time{1}, 1:nChan, out_ord);
axis tight;
colorbar;
xlabel('Time (s)');
ylabel('Channels (grouped by ROI)');
title('Channel x Time grouped by Neuromag306 ROIs');

edges   = [find(diff(roi_idx_sorted)~=0); numel(roi_idx_sorted)];
starts  = [1; edges(1:end-1)+1];
centers = floor((starts + edges)/2);

yticks(centers);
yticklabels(roi_names);
set(gca,'TickDir','out','FontSize',10);

hold on;
for e = edges(1:end-1)'
    yline(e+0.5,'k-');
end
hold off;

