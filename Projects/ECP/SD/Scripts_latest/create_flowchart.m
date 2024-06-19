function create_flowchart()
    figure;
    hold on;
    
    % Define positions
    boxWidth = 120;
    boxHeight = 40;
    verticalSpacing = 100;
    horizontalSpacing = 150;
    
    % Node positions
    positions = [
        1, 3; % Start
        1, 2; % Task/Conditions
        0, 1; % MEG
        2, 1; % fMRI
        1, 0; % Subjects
        0, -1; % Patients
        2, -1; % Control
        1, -2; % Preprocessing
        1, -3; % Source Analysis
        1, -4; % Contrasts
        1, -5; % Compute LI
        1, -6; % Atlas and ROIs
        1, -7; % Comparisons and Interpretations
];
    
    % Scale positions
    positions(:, 1) = positions(:, 1) * horizontalSpacing;
    positions(:, 2) = positions(:, 2) * verticalSpacing;
    
    % Node labels
    labels = {
        'Dataset - MEG and fMRI', 'Task/Conditions', 'MEG: Semantic Decision', 'fMRI: Tone', ...
        'Subjects', 'Patients with Temporal Lobe Epilepsy', 'Healthy Control', 'Preprocessing', ...
        'Source Analysis', 'Contrasts', 'Compute LI', 'Atlas and ROIs', ...
        'Comparisons and Interpretations', ...
    };
    
    % Draw nodes
    for i = 1:length(labels)
        drawNode(positions(i, 1), positions(i, 2), boxWidth, boxHeight, labels{i});
    end
    
    % Draw connections
    drawArrow([positions(1, :); positions(2, :)]);
    drawArrow([positions(2, :); positions(3, :)]);
    drawArrow([positions(2, :); positions(4, :)]);
    drawArrow([positions(3, :); positions(5, :)]);
    drawArrow([positions(4, :); positions(5, :)]);
    drawArrow([positions(5, :); positions(6, :)]);
    drawArrow([positions(5, :); positions(7, :)]);
    drawArrow([positions(6, :); positions(8, :)]);
    drawArrow([positions(7, :); positions(8, :)]);
    drawArrow([positions(8, :); positions(9, :)]);
    drawArrow([positions(9, :); positions(10, :)]);
    drawArrow([positions(10, :); positions(11, :)]);
    drawArrow([positions(11, :); positions(12, :)]);
    drawArrow([positions(12, :); positions(13, :)]);
    drawArrow([positions(13, :); positions(14, :)]);
    
    hold off;
    axis equal;
    axis off;
end

function drawNode(x, y, width, height, label)
    rectangle('Position', [x-width/2, y-height/2, width, height], 'EdgeColor', 'b');
    text(x, y, label, 'HorizontalAlignment', 'center');
end

function drawArrow(pos)
    line(pos(:, 1), pos(:, 2), 'Color', 'k', 'LineWidth', 2, 'Marker', '>', 'MarkerIndices', [2]);
end
