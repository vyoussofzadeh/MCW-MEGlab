function varargout = process_export_screenshot( varargin )
% PROCESS_SNAPSHOT: Save snapshot.
%
% USAGE:     sProcess = process_snapshot('GetDescription')
%                       process_snapshot('Run', sProcess, sInputs)

% @=============================================================================
% This function is part of the Brainstorm software:
% https://neuroimage.usc.edu/brainstorm
%
% Copyright (c) University of Southern California & McGill University
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPLv3
% license can be found at http://www.gnu.org/copyleft/gpl.html.
%
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
% For more information type "brainstorm license" at command prompt.
% =============================================================================@
%
% Authors: Francois Tadel, 2012-2022; Vahab Youssof Zadeh, 2024

eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
% Description the process
sProcess.Comment     = 'Save snapshot (MCW_clinical)';
sProcess.Category    = 'Custom';
sProcess.SubGroup    = 'File';
sProcess.Index       = 982;
sProcess.Description = 'https://neuroimage.usc.edu/brainstorm/Tutorials/Scripting';
% Definition of the input accepted by this process
sProcess.InputTypes  = {'raw', 'data', 'results', 'timefreq', 'matrix', 'dipoles', 'pdata', 'presults', 'ptimefreq', 'pmatrix'};
sProcess.OutputTypes = {'raw', 'data', 'results', 'timefreq', 'matrix', 'dipoles', 'pdata', 'presults', 'ptimefreq', 'pmatrix'};
sProcess.nInputs     = 1;
sProcess.nMinFiles   = 1;
% === Orientation Options
sProcess.options.orientpreset.Comment = 'Orientation Preset: ';
sProcess.options.orientpreset.Type    = 'combobox';
sProcess.options.orientpreset.Value   = {1, {
    'left;right;top;bottom;left_intern;right_intern', ...
    'bottom;top;left;right;left_intern;right_intern', ...
    'left;right;top;bottom;right_intern;left_intern', ...
    'left;bottom;right', ...
    'left;right;top;bottom', ...
    'Custom (Specify Below)'
    }, 'Type', 'Label', 'Label', 'Value'};

% === Custom Orientation
% sProcess.options.customorient.Comment = 'Custom Orientation (e.g., {"left", "right", "top"}): ';
sProcess.options.customorient.Comment = 'Custom Orientation (e.g., left, right, top): ';
sProcess.options.customorient.Type    = 'text';
sProcess.options.customorient.Value   = '';

% === Background Color Selection
sProcess.options.background.Comment = 'Background color: ';
sProcess.options.background.Type    = 'combobox';
sProcess.options.background.Value   = {1, {'White', 'Black'}, 'Type', 'Label', 'Label', 'Value'};

% === TIME: Single view
sProcess.options.time.Comment = 'Time (in seconds):';
sProcess.options.time.Type    = 'value';
sProcess.options.time.Value   = {0, 's', 4};
% === TIME: Contact sheet
sProcess.options.contact_time.Comment = 'Contact sheet (start time, stop time):';
sProcess.options.contact_time.Type    = 'value';
sProcess.options.contact_time.Value   = {[0,.1], 'list', 4};
sProcess.options.contact_nimage.Comment = 'Contact sheet (number of images):';
sProcess.options.contact_nimage.Type    = 'value';
sProcess.options.contact_nimage.Value   = {12, '', 0};
% === THRESOLD
sProcess.options.threshold.Comment = 'Amplitude threshold:';
sProcess.options.threshold.Type    = 'value';
sProcess.options.threshold.Value   = {30, '%', 0};
% === ROW NAMES
sProcess.options.rowname.Comment    = 'Time-frequency signal name (empty=all): ';
sProcess.options.rowname.Type       = 'text';
sProcess.options.rowname.Value      = '';
sProcess.options.rowname.InputTypes = {'timefreq', 'matrix'};

% === COMMENT
% Add a folder directory input
sProcess.options.savedir.Comment = 'Saving Dir.:';
sProcess.options.savedir.Type    = 'text';
sProcess.options.savedir.Value   = ''; % Default value can be empty or a specific path

% Add a saving file name input
sProcess.options.sname.Comment = 'Saving filename:';
sProcess.options.sname.Type    = 'text';
sProcess.options.sname.Value   = ''; % Default value can be empty or a specific name

end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
iType = strcmpi(sProcess.options.type.Value{1}, sProcess.options.type.Value{2}(2,:));
Comment = ['Snapshot: ' sProcess.options.type.Value{2}{1,iType}];
end

%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
% Returned files: same as input
OutputFiles = {sInputs.FileName};
% Get options

% Get Orientation Selection
orientPreset = sProcess.options.orientpreset.Value{1};
customOrient = sProcess.options.customorient.Value;

% Determine Orientations based on selection
switch orientPreset
    case 1
        Orient = {'left', 'right', 'top', 'bottom', 'left_intern', 'right_intern'};
    case 2
        Orient = {'bottom', 'top', 'left', 'right', 'left_intern', 'right_intern'};
    case 3
        Orient = {'left', 'right', 'top', 'bottom', 'right_intern', 'left_intern'};
    case 4
        Orient = {'left', 'bottom', 'right'};
    case 5
        Orient = {'left', 'right', 'top', 'bottom'};
    case 6
        % Custom Orientation
        if ~isempty(customOrient)
            Orient = strsplit(customOrient, ', ');  % Splits at the comma and space
        else
            bst_report('Error', sProcess, [], 'Custom orientation is selected but not defined.');
            return;
        end
end

%%
% Get Background Color Selection
backgroundSelection = sProcess.options.background.Value{1};

fname = OutputFiles{1};

% Obtain saving directory
savedir = sProcess.options.savedir.Value;
svname = sProcess.options.sname.Value;

close all
hFig = view_surface_data([], fname, [], 'NewFigure');

% Set background color based on user selection
if backgroundSelection == 1
    set(hFig, 'color', 'w');
elseif backgroundSelection == 2
    set(hFig, 'color', 'k');
else
    % Default or error handling
    disp('Invalid background selection. Using default white background.');
    set(hFig, 'color', 'w');
end

bst_colormaps('SetColorbarVisible', hFig, 0);
axis equal
pause,

%% Export images: PNG
b = []; a = [];
img = [];
for iOrient=1:length(Orient)
    figure_3d('SetStandardView', hFig, Orient{iOrient});
    img{iOrient} = out_figure_image(hFig, '', '');
    b=double(img{iOrient});
    if isa(b,'uint8'), b=double(b)/255; end
    if max(b(:))>1, b=double(b)/double(max(b(:))); end
    a{iOrient}=double(b);
    pause(1)
end

ncut = 16;
domosaiccrop = 1;
if domosaiccrop
    cropt_idx={};
    for n=1:numel(a)
        cropt=any(any(diff(a{n},1,2),2),3);
        cropt_idx{n,1}=max(1,sum(~cumsum(cropt))-ncut):size(a{n},1)-max(0,sum(~cumsum(flipud(cropt)))-ncut);
        cropt=any(any(diff(a{n},1,1),1),3);
        cropt_idx{n,2}=max(1,sum(~cumsum(cropt))-ncut):size(a{n},2)-max(0,sum(~cumsum(flipud(cropt)))-ncut);
    end
end

if domosaiccrop
    cropt_idx1234=cropt_idx{1,1}; for n=2:numel(a), cropt_idx1234=union(cropt_idx1234,cropt_idx{n,1}); end
    ta=[]; for n=1:numel(a), ta=cat(2,ta,a{n}(cropt_idx1234,cropt_idx{n,2},:)); end; a=ta;
else
    a=cat(2,a{:});
end

imgFile = fullfile(savedir, [svname,'.png']);
imwrite(a,imgFile);

% -PNG
combined_path = fullfile(savedir,[svname, '.png']); 
web(combined_path, '-new');

%% -SVG (incomplete - adding extra white border)
% Create a figure
% fig = figure;
% imshow(a);  % Display the image
% set(gca, 'Position', [0 0 1 1]); % Optional: Remove any margins
% set(gca, 'color', 'none');
% combined_path = fullfile(savedir, [svname, '.svg']);
% print(fig, '-dsvg', combined_path);  % Save the figure as SVG
% 
% % Close the figure
% close(fig);

end