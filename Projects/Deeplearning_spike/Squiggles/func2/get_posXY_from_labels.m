function posXY = get_posXY_from_labels(refLbl)
% Returns [C x 2] XY positions in [-1,1] for Neuromag 306 labels.

% 1) Prep FieldTrip layout
cfg = [];
cfg.layout = 'neuromag306all.lay';          % ships with FieldTrip
lay = ft_prepare_layout(cfg);

% 2) Map label -> position
lab2idx = containers.Map(lay.label, 1:numel(lay.label));
C = numel(refLbl);
posXY = nan(C,2);
for i = 1:C
    li = refLbl{i};
    if isKey(lab2idx, li)
        j = lab2idx(li);
        posXY(i,:) = lay.pos(j,1:2);        % [x y]
    else
        % some datasets use 'MEG011' instead of 'MEG0111'; try best-effort match
        k = find(strcmp(lay.label, li(1:min(end,8))), 1);
        if ~isempty(k)
            posXY(i,:) = lay.pos(k,1:2);
        else
            warning('No layout position for %s; setting [0 0].', li);
            posXY(i,:) = [0 0];
        end
    end
end

% 3) Normalize to [-1,1]
posXY = normalize(posXY, "range", [-1 1]);
end
