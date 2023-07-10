function varargout = process_merge_dipoles( varargin )
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
sProcess.Comment     = 'Merge BS dipoles';
sProcess.Category    = 'Custom';
sProcess.SubGroup    = 'Sources';
sProcess.Index       = 334;
sProcess.Description = 'https://neuroimage.usc.edu/brainstorm/Tutorials/CoregisterSubjects';
% Definition of the input accepted by this process
% sProcess.InputTypes  = {'results','data', 'raw'};
sProcess.InputTypes  = {'results','dipoles'};
sProcess.OutputTypes = {'results'};
% sProcess.InputTypes  = {'data', 'results', 'timefreq', 'matrix'};
% sProcess.OutputTypes = {'data', 'results', 'timefreq', 'matrix'};
sProcess.nInputs     = 1;
% sProcess.nMinFiles   = 1;
sProcess.nMinFiles   = 2;

end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
Comment = sProcess.Comment;
end

%% ===== RUN =====
function OutputFiles = Run(sProcess, sInput) %#ok<DEFNU>
OutputFiles = {};

%%
cd_org = pwd;
k=1;
item = [];
D=[];
for i=1:length(sInput)
    sResultP = in_bst_results(sInput(i).FileName, 0);
    if contains(sResultP.Comment, 'Dipoles')
        ProtocolInfo =        bst_get('ProtocolInfo');
        DipolesFile = bst_get('DipolesForFile',    sInput(i).FileName, sInput(i).iStudy);
        [pathstr, ~, ~] = fileparts(DipolesFile.FileName);
        cd(fullfile(ProtocolInfo.STUDIES, pathstr))
        if isempty(D)
            D = rdir('dipoles_fit*.mat');
        end
        item(k) = sInput(i).iItem;
        k =k+1;
    end
end

DipoleFiles = [];
if ~isempty(D)
    for j=1:length(D)
        DipoleFiles{1,j} = fullfile(D(j).folder, D(j).name);
    end
end

DipoleFiles = DipoleFiles(item);
disp(DipoleFiles')

cd(cd_org)

%%
% Load in all dipole files
for i = 1:length(DipoleFiles)
    % Check file path file
    if ~file_exist(DipoleFiles{i})
        DipoleFiles{i} = file_fullpath(DipoleFiles{i});
    end
    if ~file_exist(DipoleFiles{i})
        error(['File not found: ' DipoleFiles{i}]);
    end
    % Load file
    DipoleMat = load(DipoleFiles{i});
    % Do not accept subsets of dipoles
    if (DipoleMat.Subset > 1)
        error('TODO: Update this function for merging files with subsets.')
    end
    % First file: template
    if (i == 1)
        MergeMat = DipoleMat;
        MergeMat.DipoleNames = {};
        MergeMat.Dipole = [];
        
        % Calculate sampling rate of first dipole file
        if (length(DipoleMat.Time) == 1)
            DipOneSamplingRate = 0;  % Set sampling rate to 0 if only one dipole
        elseif (length(DipoleMat.Time) > 1)
            DipOneSamplingRate = DipoleMat.Time(2) - DipoleMat.Time(1);
        end
        
        % Following files: check that they are compatible
    else
        if ((DipOneSamplingRate ~= 0) && (length(DipoleMat.Time) > 1))
            if (abs(DipOneSamplingRate - (DipoleMat.Time(2) - DipoleMat.Time(1))) > 0.00001)
                error('Only files with equal sampling rate can be merged.');
            end
        elseif (DipoleMat.Subset ~= MergeMat.Subset)
            error('Only files with equal number of channel subsets can be merged.');
        end
    end
    
    % Names just have the group number (i.e. Group #1, Group #2)
    lastGroupNumber = length(MergeMat.DipoleNames);
    nGroups = length(DipoleMat.DipoleNames);
    % Loop through the groups
    for g = 1:nGroups
        % Add the name
        MergeMat.DipoleNames{end+1} = ['Group #' num2str(lastGroupNumber + g)];
        % Update the index numbers
        ind = find([DipoleMat.Dipole.Index] == g);
        if ~isempty(ind)
            [DipoleMat.Dipole(ind).Index] = deal(lastGroupNumber + g);
        end
    end
    % Merge the Dipole structures
    MergeMat.Dipole = [MergeMat.Dipole, DipoleMat.Dipole];
    MergeMat.DataFile = '';
    % Merge history
    MergeMat = bst_history('add', MergeMat, 'merge', ['Merged file: ' file_short(DipoleFiles{i})]);
end

% Set new Time sampling for merged dipoles
MergeMat.Time = unique([MergeMat.Dipole.Time]);

% ===== SAVE NEW FILE =====
% Update the comment to reflect the number of merged files
MergeMat.Comment = ['Merge: ' num2str(length(DipoleFiles)) ' files'];
% Create output filename
OutputFile = file_unique(bst_fullfile(fileparts(DipoleFiles{1}), 'dipoles_merged.mat'));
% Save new file in Brainstorm format
bst_save(OutputFile, MergeMat, 'v7');


% ===== UPDATE DATABASE =====
% Get study of the first file
[sStudy,iStudy] = bst_get('DipolesFile', file_short(DipoleFiles{1}));
% Create structure
BstDipolesMat = db_template('Dipoles');
BstDipolesMat.FileName = file_short(OutputFile);
BstDipolesMat.Comment  = MergeMat.Comment;
BstDipolesMat.DataFile = MergeMat.DataFile;
% Add to study
sStudy = bst_get('Study', iStudy);
iDipole = length(sStudy.Dipoles) + 1;
sStudy.Dipoles(iDipole) = BstDipolesMat;
% Save study
bst_set('Study', iStudy, sStudy);
% Update tree
panel_protocols('UpdateNode', 'Study', iStudy);
% Select node
panel_protocols('SelectNode', [], BstDipolesMat.FileName);
% Save database
db_save();

end

function [varargout] = rdir(rootdir,varargin)

% use the current directory if nothing is specified
if ~exist('rootdir','var')
    rootdir = '*';
end

% split the file path around the wild card specifiers
prepath = '';       % the path before the wild card
wildpath = '';      % the path wild card
postpath = rootdir; % the path after the wild card
I = find(rootdir==filesep,1,'last');
if ~isempty(I)
    prepath = rootdir(1:I);
    postpath = rootdir(I+1:end);
    I = find(prepath=='*',1,'first');
    if ~isempty(I)
        postpath = [prepath(I:end) postpath];
        prepath = prepath(1:I-1);
        I = find(prepath==filesep,1,'last');
        if ~isempty(I)
            wildpath = prepath(I+1:end);
            prepath = prepath(1:I);
        end
        I = find(postpath==filesep,1,'first');
        if ~isempty(I)
            wildpath = [wildpath postpath(1:I-1)];
            postpath = postpath(I:end);
        end
    end
end
% disp([' "' prepath '" ~ "' wildpath '" ~ "' postpath '" ']);


if isempty(wildpath)
    % if no directory wildcards then just get file list
    D = dir([prepath postpath]);
    D([D.isdir]==1) = [];
    for ii = 1:length(D)
        if (~D(ii).isdir)
            D(ii).name = [prepath D(ii).name];
        end
    end
    
    % disp(sprintf('Scanning "%s"   %g files found',[prepath postpath],length(D)));
    
elseif strcmp(wildpath,'**') % a double wild directory means recurs down into sub directories
    
    % first look for files in the current directory (remove extra filesep)
    D = rdir([prepath postpath(2:end)]);
    
    % then look for sub directories
    Dt = dir('');
    tmp = dir([prepath '*']);
    % process each directory
    for ii = 1:length(tmp)
        if (tmp(ii).isdir && ~strcmpi(tmp(ii).name,'.') && ~strcmpi(tmp(ii).name,'..') ),
            Dt = [Dt; rdir([prepath tmp(ii).name filesep wildpath postpath])];
        end
    end
    D = [D; Dt];
    
else
    % Process directory wild card looking for sub directories that match
    tmp = dir([prepath wildpath]);
    D = dir('');
    % process each directory found
    for ii = 1:length(tmp)
        if (tmp(ii).isdir && ~strcmpi(tmp(ii).name,'.') && ~strcmpi(tmp(ii).name,'..') ),
            D = [D; rdir([prepath tmp(ii).name postpath])];
        end
    end
end


% Apply filter
if (nargin>=2 && ~isempty(varargin{1}))
    date = [D.date];
    datenum = [D.datenum];
    bytes = [D.bytes];
    
    try
        eval(sprintf('D((%s)==0) = [];',varargin{1}));
    catch
        warning('Error: Invalid TEST "%s"',varargin{1});
    end
end

% display listing if no output variables are specified
if nargout==0
    pp = {'' 'k' 'M' 'G' 'T'};
    for ii=1:length(D)
        sz = D(ii).bytes;
        if sz<=0
            disp(sprintf(' %31s %-64s','',D(ii).name));
        else
            ss = min(4,floor(log2(sz)/10));
            disp(sprintf('%4.0f %1sb   %20s   %-64s ',sz/1024^ss,pp{ss+1},D(ii).date,D(ii).name));
        end
    end
else
    % send list out
    varargout{1} = D;
end

end

function varargout = process_merge_dipoles( varargin )
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
sProcess.Comment     = 'Merge BS dipoles';
sProcess.Category    = 'Custom';
sProcess.SubGroup    = 'Sources';
sProcess.Index       = 334;
sProcess.Description = 'https://neuroimage.usc.edu/brainstorm/Tutorials/CoregisterSubjects';
% Definition of the input accepted by this process
% sProcess.InputTypes  = {'results','data', 'raw'};
sProcess.InputTypes  = {'results','dipoles'};
sProcess.OutputTypes = {'results'};
% sProcess.InputTypes  = {'data', 'results', 'timefreq', 'matrix'};
% sProcess.OutputTypes = {'data', 'results', 'timefreq', 'matrix'};
sProcess.nInputs     = 1;
% sProcess.nMinFiles   = 1;
sProcess.nMinFiles   = 2;

end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
Comment = sProcess.Comment;
end

%% ===== RUN =====
function OutputFiles = Run(sProcess, sInput) %#ok<DEFNU>
OutputFiles = {};

%%
cd_org = pwd;
k=1;
item = [];
D=[];
for i=1:length(sInput)
    sResultP = in_bst_results(sInput(i).FileName, 0);
    if contains(sResultP.Comment, 'Dipoles')
        ProtocolInfo =        bst_get('ProtocolInfo');
        DipolesFile = bst_get('DipolesForFile',    sInput(i).FileName, sInput(i).iStudy);
        [pathstr, ~, ~] = fileparts(DipolesFile.FileName);
        cd(fullfile(ProtocolInfo.STUDIES, pathstr))
        if isempty(D)
            D = rdir('dipoles_fit*.mat');
        end
        item(k) = sInput(i).iItem;
        k =k+1;
        %         break,
    end
end

DipoleFiles = [];
if ~isempty(D)
    for j=1:length(D)
        DipoleFiles{1,j} = fullfile(D(j).folder, D(j).name);
    end
end

DipoleFiles = DipoleFiles(item);
disp(DipoleFiles')

cd(cd_org)

%%
% Load in all dipole files
for i = 1:length(DipoleFiles)
    % Check file path file
    if ~file_exist(DipoleFiles{i})
        DipoleFiles{i} = file_fullpath(DipoleFiles{i});
    end
    if ~file_exist(DipoleFiles{i})
        error(['File not found: ' DipoleFiles{i}]);
    end
    % Load file
    DipoleMat = load(DipoleFiles{i});
    % Do not accept subsets of dipoles
    if (DipoleMat.Subset > 1)
        error('TODO: Update this function for merging files with subsets.')
    end
    % First file: template
    if (i == 1)
        MergeMat = DipoleMat;
        MergeMat.DipoleNames = {};
        MergeMat.Dipole = [];
        
        % Calculate sampling rate of first dipole file
        if (length(DipoleMat.Time) == 1)
            DipOneSamplingRate = 0;  % Set sampling rate to 0 if only one dipole
        elseif (length(DipoleMat.Time) > 1)
            DipOneSamplingRate = DipoleMat.Time(2) - DipoleMat.Time(1);
        end
        
        % Following files: check that they are compatible
    else
        if ((DipOneSamplingRate ~= 0) && (length(DipoleMat.Time) > 1))
            if (abs(DipOneSamplingRate - (DipoleMat.Time(2) - DipoleMat.Time(1))) > 0.00001)
                error('Only files with equal sampling rate can be merged.');
            end
        elseif (DipoleMat.Subset ~= MergeMat.Subset)
            error('Only files with equal number of channel subsets can be merged.');
        end
    end
    
    % Names just have the group number (i.e. Group #1, Group #2)
    lastGroupNumber = length(MergeMat.DipoleNames);
    nGroups = length(DipoleMat.DipoleNames);
    % Loop through the groups
    for g = 1:nGroups
        % Add the name
        MergeMat.DipoleNames{end+1} = ['Group #' num2str(lastGroupNumber + g)];
        % Update the index numbers
        ind = find([DipoleMat.Dipole.Index] == g);
        if ~isempty(ind)
            [DipoleMat.Dipole(ind).Index] = deal(lastGroupNumber + g);
        end
    end
    % Merge the Dipole structures
    MergeMat.Dipole = [MergeMat.Dipole, DipoleMat.Dipole];
    MergeMat.DataFile = '';
    % Merge history
    MergeMat = bst_history('add', MergeMat, 'merge', ['Merged file: ' file_short(DipoleFiles{i})]);
end

% Set new Time sampling for merged dipoles
MergeMat.Time = unique([MergeMat.Dipole.Time]);

% ===== SAVE NEW FILE =====
% Update the comment to reflect the number of merged files
MergeMat.Comment = ['Merge: ' num2str(length(DipoleFiles)) ' files'];
% Create output filename
OutputFile = file_unique(bst_fullfile(fileparts(DipoleFiles{1}), 'dipoles_merged.mat'));
% Save new file in Brainstorm format
bst_save(OutputFile, MergeMat, 'v7');


% ===== UPDATE DATABASE =====
% Get study of the first file
[sStudy,iStudy] = bst_get('DipolesFile', file_short(DipoleFiles{1}));
% Create structure
BstDipolesMat = db_template('Dipoles');
BstDipolesMat.FileName = file_short(OutputFile);
BstDipolesMat.Comment  = MergeMat.Comment;
BstDipolesMat.DataFile = MergeMat.DataFile;
% Add to study
sStudy = bst_get('Study', iStudy);
iDipole = length(sStudy.Dipoles) + 1;
sStudy.Dipoles(iDipole) = BstDipolesMat;
% Save study
bst_set('Study', iStudy, sStudy);
% Update tree
panel_protocols('UpdateNode', 'Study', iStudy);
% Select node
panel_protocols('SelectNode', [], BstDipolesMat.FileName);
% Save database
db_save();

end

function [varargout] = rdir(rootdir,varargin)

% use the current directory if nothing is specified
if ~exist('rootdir','var')
    rootdir = '*';
end

% split the file path around the wild card specifiers
prepath = '';       % the path before the wild card
wildpath = '';      % the path wild card
postpath = rootdir; % the path after the wild card
I = find(rootdir==filesep,1,'last');
if ~isempty(I)
    prepath = rootdir(1:I);
    postpath = rootdir(I+1:end);
    I = find(prepath=='*',1,'first');
    if ~isempty(I)
        postpath = [prepath(I:end) postpath];
        prepath = prepath(1:I-1);
        I = find(prepath==filesep,1,'last');
        if ~isempty(I)
            wildpath = prepath(I+1:end);
            prepath = prepath(1:I);
        end
        I = find(postpath==filesep,1,'first');
        if ~isempty(I)
            wildpath = [wildpath postpath(1:I-1)];
            postpath = postpath(I:end);
        end
    end
end
% disp([' "' prepath '" ~ "' wildpath '" ~ "' postpath '" ']);


if isempty(wildpath)
    % if no directory wildcards then just get file list
    D = dir([prepath postpath]);
    D([D.isdir]==1) = [];
    for ii = 1:length(D)
        if (~D(ii).isdir)
            D(ii).name = [prepath D(ii).name];
        end
    end
    
    % disp(sprintf('Scanning "%s"   %g files found',[prepath postpath],length(D)));
    
elseif strcmp(wildpath,'**') % a double wild directory means recurs down into sub directories
    
    % first look for files in the current directory (remove extra filesep)
    D = rdir([prepath postpath(2:end)]);
    
    % then look for sub directories
    Dt = dir('');
    tmp = dir([prepath '*']);
    % process each directory
    for ii = 1:length(tmp)
        if (tmp(ii).isdir && ~strcmpi(tmp(ii).name,'.') && ~strcmpi(tmp(ii).name,'..') ),
            Dt = [Dt; rdir([prepath tmp(ii).name filesep wildpath postpath])];
        end
    end
    D = [D; Dt];
    
else
    % Process directory wild card looking for sub directories that match
    tmp = dir([prepath wildpath]);
    D = dir('');
    % process each directory found
    for ii = 1:length(tmp)
        if (tmp(ii).isdir && ~strcmpi(tmp(ii).name,'.') && ~strcmpi(tmp(ii).name,'..') ),
            D = [D; rdir([prepath tmp(ii).name postpath])];
        end
    end
end


% Apply filter
if (nargin>=2 && ~isempty(varargin{1}))
    date = [D.date];
    datenum = [D.datenum];
    bytes = [D.bytes];
    
    try
        eval(sprintf('D((%s)==0) = [];',varargin{1}));
    catch
        warning('Error: Invalid TEST "%s"',varargin{1});
    end
end

% display listing if no output variables are specified
if nargout==0
    pp = {'' 'k' 'M' 'G' 'T'};
    for ii=1:length(D)
        sz = D(ii).bytes;
        if sz<=0
            disp(sprintf(' %31s %-64s','',D(ii).name));
        else
            ss = min(4,floor(log2(sz)/10));
            disp(sprintf('%4.0f %1sb   %20s   %-64s ',sz/1024^ss,pp{ss+1},D(ii).date,D(ii).name));
        end
    end
else
    % send list out
    varargout{1} = D;
end

end