function varargout = process_ROI_Peak_Coordinate(varargin )
% PROCESS_FT_SOURCEANALYSIS Call FieldTrip function ft_sourceanalysis (DICS)

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
% Authors: Vahab YoussofZadeh, 2023

eval(macro_method);
end

%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>

% Description of the process
sProcess.Comment     = 'ROI peak coordinate';
sProcess.Category    = 'Custom';
sProcess.SubGroup    = 'Sources';
sProcess.Index       = 337;
sProcess.Description = 'https://neuroimage.usc.edu/brainstorm/Tutorials/CoregisterSubjects';
sProcess.InputTypes  = {'results'};
sProcess.OutputTypes = {'results'};
sProcess.nInputs     = 1;
sProcess.nMinFiles   = 1;

% Add a number of peaks input
sProcess.options.numpeaks.Comment = 'Number of peaks:';
sProcess.options.numpeaks.Type    = 'value';
sProcess.options.numpeaks.Value   = {1, '', 0}; % Default value is 1, second parameter is unit (none in this case), third is minimum number of peaks

% Add a checkbox for whether to save the data
sProcess.options.saveData.Comment = 'Save data?';
sProcess.options.saveData.Type    = 'checkbox';
sProcess.options.saveData.Value   = 1; % Default checked

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
Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInput)

OutputFiles = {};
sResultP = in_bst_results(sInput.FileName, 1);

% Obtain saving directory
npeaks = sProcess.options.numpeaks.Value{1};

fname = sInput.FileName;
SourceMat = in_bst(fname);

%%
sVolHeadModelFileName = in_bst_headmodel(sResultP.HeadModelFile, [], 'SurfaceFile', 'GridLoc');
sSurf = in_tess_bst(sVolHeadModelFileName.SurfaceFile);
sSubject = bst_get('SurfaceFile', sVolHeadModelFileName.SurfaceFile);
sMri = in_mri_bst(sSubject.Anatomy(sSubject.iAnatomy).FileName);
ImageGridAmp = SourceMat.ImageGridAmp;

%%
[maxval, maxidx] = max(ImageGridAmp);

b = [];
b.Vertices = sSurf.Vertices;
b.Faces = sSurf.Faces;
% figure; ft_plot_mesh(b, 'edgecolor', 'none', 'facealpha', 0.4);
% view([-180,0])
% hold on
% plot3(b.Vertices(maxidx,1),b.Vertices(maxidx,2),b.Vertices(maxidx,3),'*')

%%
coor = 1000.*b.Vertices;
[~, sor_idx] = sort(ImageGridAmp,'descend');

%%
seed_coor = coor(sor_idx(1),:);

hFig = view_surface_data([], fname, [], 'NewFigure');
set(hFig,'color','w');
bst_colormaps('SetColorbarVisible', hFig, 0);
axis equal
alpha(0.5)

seed_val = [];
marker_size = 80;
hold on
x = b.Vertices(sor_idx(1),1); y = b.Vertices(sor_idx(1),2); z = b.Vertices(sor_idx(1),3);
kk = 1;
clear fh
fh(kk) = plot3(x,y,z,'m.','MarkerSize',marker_size);

r1 = x;
str = sprintf('%s', num2str(size(seed_coor,1)));

h(kk)   = text(x, y, z, str, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Interpreter', 'none','FontSize',10);
seed_val(kk) = ImageGridAmp(sor_idx(1));

for i=1:length(ImageGridAmp)-1
    if size(seed_coor,1) < npeaks
        d2 = [];
        for j=1:size(seed_coor,1)
            d2(j) = pdist2(seed_coor(j,:),coor(sor_idx(i+1),:));
        end
        md2 = min(d2);
        if md2 > 30
            seed_coor = [seed_coor;coor(sor_idx(i+1),:)];
            x = b.Vertices(sor_idx(i+1),1); y = b.Vertices(sor_idx(i+1),2); z = b.Vertices(sor_idx(i+1),3);
            kk = kk+1;
            seed_val(kk) = ImageGridAmp(sor_idx(i+1));
            fh(kk) = plot3(x,y,z,'m.','MarkerSize',marker_size);
            str = sprintf('%s', num2str(size(seed_coor,1)));
            h(kk)   = text(x, y, z, str, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Interpreter', 'none', 'FontSize',10);
        end
    end
end

%% Set anatomy
sTemplates = bst_get('AnatomyDefaults');
iTemplate = strcmpi('ICBM152', {sTemplates.Name});
sMRI = load(fullfile(sTemplates(iTemplate).FilePath,'subjectimage_T1.mat'));

%% Coordniates
P_mni = [];
for k=1:length(seed_coor)
    P_mni(k,:) = 1e3*cs_convert(sMRI, 'scs', 'mni', seed_coor(k,:)/1e3);
end

% disp('mni coordinates')
% disp(round(P_mni))
rP_mni = round(P_mni);
% disp('value')
% disp(seed_val');

%%
% Assuming P_mni and seed_val are already defined
% Round MNI coordinates
roundedMni = round(P_mni);

% Round seed values to two decimal places
roundedSeedVal = round(seed_val'*100) / 100;

% Create an index column (1, 2, 3, ...)
numValues = length(roundedMni); % Number of rows in the MNI data
index = (1:numValues)'; % Column vector of indices

% Create the table with named columns and the index
T = table(index, roundedSeedVal, roundedMni, 'VariableNames', {'index', 'value', 'mni'});

display(T)

% Check if the user wants to save the data
if sProcess.options.saveData.Value
    savedir = sProcess.options.savedir.Value;
    svname = sProcess.options.sname.Value;
    
    disp(fullfile(savedir))
    disp('as,')
    disp([svname,'.tif'])
    
    writetable(T,[fullfile(savedir,svname),'.txt']);
    
end
%% left and right coordinates,
idx_left = find(rP_mni(:,1)<0);
idx_right = find(rP_mni(:,1)>0);

%% update color,
for j=1:length(idx_left), set(fh(idx_left(j)),'Color',[.98 .45 .02]),end
for j=1:length(idx_right), set(fh(idx_right(j)),'Color',[.88 .15 .42]),end

%% Apply trasnparency
% Set marker face color for left hemisphere
for j = 1:length(idx_left)
    set(fh(idx_left(j)), 'MarkerFaceColor', 'r'); % Set to desired color
end

% Set marker face color for right hemisphere
for j = 1:length(idx_right)
    set(fh(idx_right(j)), 'MarkerFaceColor', 'g'); % Set to desired color
end


%%
bst_progress('stop');

end