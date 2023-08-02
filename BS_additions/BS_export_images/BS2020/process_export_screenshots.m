function varargout = process_export_screenshots( varargin )
% PROCESS_DIPOLE_SCANNING: Generates a brainstorm dipole file from the GLS and GLS-P inverse solutions.

% @=============================================================================
% This function is part of the Brainstorm software:
% https://neuroimage.usc.edu/brainstorm
%
% Copyright (c)2000-2020 University of Southern California & McGill University
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
% Authors: Vahab Youssof Zadeh, 2023

eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>

% Description the process
sProcess.Comment     = 'Export screenshot (MEG reports)';
sProcess.Category    = 'Custom';
sProcess.SubGroup    = 'File';
sProcess.Index       = 981;
sProcess.Description = 'https://neuroimage.usc.edu/brainstorm/Tutorials/CoregisterSubjects';
sProcess.InputTypes  = {'results'};
sProcess.OutputTypes = {'results'};
sProcess.nInputs     = 1;
sProcess.nMinFiles   = 1;

end

%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
Comment = sProcess.Comment;
end

%% ===== RUN =====
function OutputFiles = Run(~, sInput) %#ok<DEFNU>
OutputFiles = {};

sResultP = in_bst_results(sInput.FileName, 1);
ProtocolInfo =        bst_get('ProtocolInfo');
[sStudy, ~, ~] = bst_get('AnyFile', sInput.FileName);


FileName = sInput.FileName;

if isempty(sStudy)
    [~, fileBase, fileExt] = bst_fileparts(FileName);
    fileBase = [fileBase, fileExt];
else
    ProtocolInfo = bst_get('ProtocolInfo');
    [FileName, ~, isAnatomy] = file_fullpath(FileName);
    if isAnatomy
        filePath = ProtocolInfo.SUBJECTS;
    else
        filePath = ProtocolInfo.STUDIES;
    end
    fileBase = file_win2unix(strrep(FileName, filePath, ''));
end

%%
disp('Where the map is goingn to be exported to ...')
svdir = input('1: enter saving dir:','s');svdir = strrep(svdir, ' ', '');

BSpath = ProtocolInfo.STUDIES;
fname = fileBase;
cd(BSpath)

%%
sfile = load(fname);
svname = sfile.Comment;
disp(['suggesting name:',svname]);
svname = input('enter saving name:', 's');
disp('===========')


clc
disp('1: left;right;top;bottom;left_intern;right_intern')
disp('2: bottom;top;left;right;left_intern;right_intern')
disp('3: left;right;top;bottom;right_intern;left_intern')
disp('4: left;bottom;right')
disp('5: left;right;top;bottom')
disp('6: optional, e.g, {left;right;top}')
disp('select views')
side_sel = input(':');

switch side_sel
    case 1
        Orient = {'left'; 'right';'top';'bottom';'left_intern';'right_intern'};
    case 2
        Orient = {'bottom';'top';'left';'right';'left_intern';'right_intern'};
    case 3
        Orient = {'left'; 'right';'top';'bottom';'right_intern';'left_intern'};
    case 4
        Orient = {'left'; 'bottom';'right'};
    case 5
        Orient = {'left'; 'right';'top';'bottom'};
    case 6
        disp('enter selected views, in quotation marks')
        side_sel_man = input('');
        Orient = side_sel_man;
end


clc
disp('1: White background')
disp('2: Black background')
backg_sel = input(':');

disp('adjust the surface theshold (in BS GUI)')
disp('then hit enter')

close all
hFig = view_surface_data([], fname, [], 'NewFigure');
switch backg_sel
    case 1
        set(hFig,'color','w');
end
bst_colormaps('SetColorbarVisible', hFig, 0);
axis equal
pause,

a = [];
for i=1:length(Orient)
    figure_3d('SetStandardView', hFig, Orient{i});
    img = out_figure_image(hFig, '', '');
    imgFile = fullfile(svdir, [svname,'.tif']);
    out_image(imgFile, img);
    b=imread(imgFile);
    if isa(b,'uint8'), b=double(b)/255; end
    if max(b(:))>1, b=double(b)/double(max(b(:))); end
    a{i}=double(b);
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

imwrite(a,imgFile);
disp('Completed!, images were saved at,')
disp(fullfile(svdir, svdir))
cd(svdir)

bst_progress('stop');

end

