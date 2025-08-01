function doPlotExport(plot_option, save_dir, figName, fileType)
% doPlotExport Conditionally exports and displays a figure
%
% USAGE:
%   doPlotExport(plot_option, save_dir, figName)
%   doPlotExport(plot_option, save_dir, figName, fileType)
%
% INPUTS:
%   plot_option : If 1, the figure is exported; otherwise, no action.
%   save_dir    : Directory where the exported figure is saved.
%   figName     : Filename (without extension) for the exported figure.
%   fileType    : (Optional) Format/extension, e.g., 'svg' (default 'svg').
%
% Example:
%   doPlotExport(1, 'C:\myFolder', 'AED_Count', 'png');

    if nargin < 4
        fileType = 'svg'; % Default file type
    end

    % Only proceed if plot_option == 1
    if plot_option == 1
        % 1) Build config structure
        cfg = [];
        cfg.outdir   = save_dir;
        cfg.filename = figName;
        cfg.type     = fileType;

        % 2) Export the current figure
        do_export_fig(cfg);

        % 3) Close all figures (optional)
        close all;

        % 4) Construct the full path and open the file in the default browser
        combined_path = fullfile(save_dir, [figName, '.', fileType]);
        web(combined_path, '-new');
    end
end
