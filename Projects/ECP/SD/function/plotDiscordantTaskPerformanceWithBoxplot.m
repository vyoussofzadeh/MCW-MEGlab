function plotDiscordantTaskPerformanceWithBoxplot(meanAccBySubject_Animal, meanAccBySubject_Falsefont, discordant_indices_anim, discordant_indices_symb, save_dir, LI_method_label, LI_method)
    % plotDiscordantTaskPerformanceWithBoxplot: Plot task performance for discordant subjects with boxplots
    %
    % Inputs:
    %   - meanAccBySubject_Animal: Table containing accuracy data for the animal task
    %   - meanAccBySubject_Falsefont: Table containing accuracy data for the symbol task
    %   - discordant_indices_anim: Indices of discordant subjects for the animal task
    %   - discordant_indices_symb: Indices of discordant subjects for the symbol task
    %   - save_dir: Directory to save the plot
    %   - LI_method_label: Cell array of LI method labels
    %   - LI_method: Index of the current LI method
    
    % Plotting Animal Task Performance with Boxplot
    figure;
    daboxplot(meanAccBySubject_Animal.mean_Animal_ACC, 'groups', ones(1, numel(meanAccBySubject_Animal.mean_Animal_ACC)), 'outsymbol', 'kx', 'xtlabels', 'test', 'fill', 0);
    xlabel('All Subjects');
    ylabel('Accuracy (%)');
    set(gca, 'FontSize', 10);
    discordant_ACC_anim = meanAccBySubject_Animal.mean_Animal_ACC(discordant_indices_anim);
    hold on;
    hScatter = scatter(ones(size(discordant_ACC_anim)), discordant_ACC_anim, 'r', 'filled');
    hold off;
    title({'Task performance (anim)', 'of Discordant Subjects'});
    l = legend(hScatter, 'Discordant', 'Location', 'south');
    legendPos = l.Position;
    legendPos(1) = legendPos(1) + 0.02;
    l.Position = legendPos;
    set(gca, 'XTick', []);
    box off;
    set(gcf, 'Position', [1000, 100, 300, 300]);
    set(gca, 'color', 'none');
    
    cfg = []; 
    cfg.outdir = save_dir; 
    filename = ['2_TaskPerformace_anim_Dicordant_LIs_', LI_method_label{LI_method}]; 
    cfg.filename = filename; 
    cfg.type = 'fig'; 
    do_export_fig(cfg);

    % Plotting Symbol Task Performance with Boxplot
    figure;
    daboxplot(meanAccBySubject_Falsefont.mean_Falsefont_ACC, 'groups', ones(1, numel(meanAccBySubject_Falsefont.mean_Falsefont_ACC)), 'outsymbol', 'kx', 'xtlabels', 'test', 'fill', 0);
    xlabel('All Subjects');
    ylabel('Accuracy (%)');
    set(gca, 'FontSize', 10);
    discordant_ACC_symb = meanAccBySubject_Falsefont.mean_Falsefont_ACC(discordant_indices_symb);
    hold on;
    hScatter = scatter(ones(size(discordant_ACC_symb)), discordant_ACC_symb, 'r', 'filled');
    hold off;
    title({'Task performance (symb)', 'of Discordant Subjects'});
    l = legend(hScatter, 'Discordant', 'Location', 'south');
    legendPos = l.Position;
    legendPos(1) = legendPos(1) + 0.02;
    l.Position = legendPos;
    set(gca, 'XTick', []);
    box off;
    set(gcf, 'Position', [1000, 100, 300, 300]);
    set(gca, 'color', 'none');
    
    cfg = []; 
    cfg.outdir = save_dir; 
    filename = ['2_TaskPerformace_symb_Dicordant_LIs_', LI_method_label{LI_method}]; 
    cfg.filename = filename; 
    cfg.type = 'fig'; 
    do_export_fig(cfg);
end
