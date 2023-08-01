function varargout = process_computeLI( varargin )
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
sProcess.Comment     = 'Compute LI (counting-based)';
sProcess.Category    = 'Custom';
sProcess.SubGroup    = 'Sources';
sProcess.Index       = 337;
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

%%
temp = load ('/MEG_data/LanguageLI/scout_Desikan-Killiany_68_template.mat');

%%

disp('1: default timing')
disp('2: define new timing')
disp('3: averaged sources')
time_interval = input('Enter time interval: ');
switch time_interval
    case 1
        clc
        disp('1: CRM')
        disp('2: DFNM')
        disp('3: PN')
        
        datatype = input('Enter Datatype: ');
        switch datatype
            case 1
                timerange=[.300 .600];  % Make sure  these are less than 1 for mS
            case 2
                timerange=[.550 1.19]; %%CUstine - changed on June 5th, 2019 from [300:1.19] to [550:1.19].
            case 3
                timerange=[.200 .600]; %%default 200 - 600ms!
        end
        disp(['[ ', num2str(timerange(1)),',', num2str(timerange(2)), '] sec was selected for LI analysis'])
    case 2
        disp('Enter Timerange: eg. [.100 .900]: ')
        timerange=input('');
    case 3
        timerange=1;  % Make sure  these are less than 1 for mS
    otherwise
        disp('There is not timerange established - pick a datatype')
        error('')
end

%%
disp('1: Positve values');
disp('2: Negative values');
disp('3: Absoulte values');
effect = input('Effects? :');
disp('==========')

switch time_interval
    case 3
        samplerate = 1;
    otherwise
        samplerate = round(inv((sResultP.Time(end)-sResultP.Time(1))/length(sResultP.Time)))-1;
end

%%
TotROI=8;  %define the number of ROIs to be used for LI-calculation. % added JL 11212014
Ratio4Threshold = 0.5; % ratio that will be multiplied by the maximum dSPM value across all ROIs within the time window of [startwindow,endwindow], in order to get a threshold for cutting off non-significant voxels. Added JL@10/30/14.

%%
disp('enter threshold ratio,(default =0.5):')
Ratio4Threshold = input(''); % Added by VYZ, 11/22/2021
if isempty(Ratio4Threshold), Ratio4Threshold = 0.5; end

%%
disp('==========')
disp('Threshold type:')
disp('1: Global-max: all time and regions combined');
disp('2: Time-max: time of interest (toi) and all regions');
disp('3: Region-max: toi and regions of interests (rois)');
Threshtype = input(''); % Added by VYZ, 09/08/2022

%% load the dSPM data and use the absolute value for these dSPM values
ImageGridAmp = sResultP.ImageGridAmp;

%%
switch effect
    case 1
        ImageGridAmp = ImageGridAmp;
    case 2
        ImageGridAmp = -ImageGridAmp;
    case 3
        ImageGridAmp= abs(ImageGridAmp);
end
% Edited by VY 8/7/19, absolute was added.
% ImageGridAmp=abs(ImageGridAmp);  % added JL 10/30/14

%% define the window of interest based on sample rate ans pretrigger data imported into brainstorm; get some max values over various time windows
%t1 = floor((startwindow+pretrig)*samplerate/1000);%start window of interest
switch time_interval
    case 3
        GlobalMax=max(ImageGridAmp(:)); % max value over all time points. Added JL @ 10/30/14
        AllMax = max(ImageGridAmp(:)); % max value over the time window of interest
    otherwise
        t1=find(sResultP.Time >= timerange(1),1); % JS 11/11/15
        t2=find(sResultP.Time >= timerange(2),1); % JS 11/11/15
        AllMax = max(max(ImageGridAmp(:,t1:t2))); % max value over the time window of interest
        GlobalMax=max(max(ImageGridAmp)); % max value over all time points. Added JL @ 10/30/14
end

%% Convert DesikenKilliany scout to select scouts
ProtocolInfo =        bst_get('ProtocolInfo');
SurfaceFile = load(fullfile(ProtocolInfo.SUBJECTS, sResultP.SurfaceFile));

Scouts = [];
sScout = [];
for i=1:length(SurfaceFile.Atlas)
    if contains(SurfaceFile.Atlas(i).Name, {'Desikan-Killiany'})
        Scouts = SurfaceFile.Atlas(i).Scouts;
    end
end
sScout.Scouts = Scouts;

%% Define Rois  Added 10/1/15 JStout
l = length(sScout.Scouts);
if l ~= 68
    
    for i=1:l
        rois{i,:} = sScout.Scouts(i).Label;
    end
    for i=1:68
        rois_temp{i,:} = temp.Scouts(i).Label;
    end
    [C,IA] = setdiff(rois_temp,rois);
    
    %%
    warning('the number of anatomical regions are not identical to atlas regions');
    disp('replacing with zero ...')
    new_eScount = temp;
    k = 1;
    for i=1:68
        if i~=IA
            new_eScount.Scouts(i) = sScout.Scouts(k);k = k+1;
        else
            new_eScount.Scouts(i).Vertices = [];
            new_eScount.Scouts(i).Seed = [];
        end
    end
    sScout  = new_eScount;
end

%%
AngSmg   = [15,16,63,64];
Front    = [3,4,5,6,11,12,25,26,29,30,33,34,37,38,39,40,41,42,49,50,53,54,55,56,57,58];
LatFront = [5,6,11,12,37,38,39,40,41,42,55,56,57,58];
LatTemp  = [1,2,17,18,31,32,61,62,65,66,67,68];
PeriSyl  = [15,16,37,38,41,42,61,62,63,64];
Tanaka   = [37,38,41,42,61,62,63,64];
Temp     = [1,2,9,10,13,14,17,18,19,20,27,28,31,32,35,36,61,62,65,66,67,68];
Whole    = 1:68;

RoiLabels = {'AngSmg', 'Front','LatFront','LatTemp', 'PeriSyl', 'Tanaka','Temp','Whole'};
RoiIndices = {AngSmg, Front,LatFront,LatTemp, PeriSyl, Tanaka, Temp,Whole};

%% Load the scout file(s)& compute LI
s1='LI_';
Summ_LI=zeros(1,TotROI); % initialize the vector that summarizes the final LIs  % added JL 11212014
Summ_LI_Label='ROI Labels: '; % initialize the string that summarizes the ROI labels  % added JL 11212014
switch time_interval
    case {2; 1}
        figure
end
plot_ind=1;
LI_label_out={};

for ii = 1:8
    
    s2 = RoiLabels{ii};
    %Odd indices are left Rois
    Ltemp_region = [];
    Ltemp_label  = [];
    hemi_roi_num=length(RoiIndices{ii});
    curr_subregion=sScout.Scouts(RoiIndices{ii});
    
    k = 1;
    for i=1:2:hemi_roi_num
        Ltemp_region=[Ltemp_region,curr_subregion(i).Vertices];
        Ltemp_label{k} = curr_subregion(i).Label; k = k+1;
    end
    
    %Even indices are right Rois
    Rtemp_region = [];
    Rtemp_label  = [];
    k = 1;
    for i=2:2:hemi_roi_num
        Rtemp_region=[Rtemp_region,curr_subregion(i).Vertices];
        Rtemp_label{k} = curr_subregion(i).Label; k = k+1;
    end
    LHscout = Ltemp_region;
    RHscout = Rtemp_region;
    
    switch time_interval % modified by VY
        case 3
            %First parse the maps into separate space-times maps for each side
            LHvals = ImageGridAmp(LHscout);
            LH_max = max(max(LHvals));
            RHvals = ImageGridAmp(RHscout);
            RH_max = max(max(RHvals));
            ROIMax = max(LH_max,RH_max);
        otherwise
            %First parse the maps into separate space-times maps for each side
            LHvals = ImageGridAmp(LHscout,t1:t2);
            LH_max = max(max(LHvals));
            RHvals = ImageGridAmp(RHscout,t1:t2);
            RH_max = max(max(RHvals));
            ROIMax = max(LH_max,RH_max);
    end
    
    switch Threshtype %modified by vy@09/08/22
        case 1
            threshold = Ratio4Threshold*GlobalMax;  % dSPM threshold to get rid of non-significant voxels. Added JL@10/30/14
        case 2
            threshold = Ratio4Threshold*AllMax;  % dSPM threshold to get rid of non-significant voxels. Added JL@10/30/14
        case 3
            threshold = Ratio4Threshold*ROIMax; % Added by VY@09/08/22
    end
    
    ind_L = find(LHvals > threshold);
    ind_R = find(RHvals > threshold);
    
    % ROI count --- above threshold voxels only
    %JL@10/30/14    indHalfROImaxL = find(LHvals > ROIMax/2);%indHalfROImaxL = find(LHvals > ROIMax/2); % Not being used right now
    L_ROIcount = length(ind_L); L_count(ii) = L_ROIcount;
    %JL@10/30/14    indHalfROImaxR = find(RHvals > ROIMax/2); % Not being used right now
    R_ROIcount = length(ind_R); R_count(ii) = R_ROIcount;
    ROIcount=sum(L_ROIcount+R_ROIcount); % to report total significant voxels over space-time
    LI_ROIcount = 100*((L_ROIcount-R_ROIcount)/(L_ROIcount+R_ROIcount));
    Summ_LI(ii)=LI_ROIcount;  % added JL 11212014
    Summ_LI_Label=[Summ_LI_Label  sprintf('\t') s2];   % added JL 11212014
    LI_label_out=[LI_label_out, s2];
    
    % ROI average --- above threshold voxels only
    LHvals_aboveThreshold = LHvals(ind_L); % a 1-D matrix, no need for mean(mean()) later
    RHvals_aboveThreshold = RHvals(ind_R); % a 1-D matrix, no need for mean(mean()) later
    L_ROIavg=mean(LHvals_aboveThreshold);
    R_ROIavg=mean(RHvals_aboveThreshold);
    ROIavg = mean([L_ROIavg,R_ROIavg]);
    LI_ROIavg = 100*((L_ROIavg-R_ROIavg)/(L_ROIavg+R_ROIavg));
    
    % Run a loop to plot LIs based on space-time voxel count as a function of threshold
    k=0;
    Rng= threshold:0.2:AllMax; % Rng= threshold:1:AllMax; % modified JL@10/30/14
    for thrTmp = Rng
        ind = find(LHvals > thrTmp);
        L_ROIcount  = length(ind);
        ind = find(RHvals > thrTmp);
        R_ROIcount  = length(ind);
        k=k+1;
        Thrshd_LI_ROIcount(k,1)=thrTmp;
        if L_ROIcount+R_ROIcount~=0
            LI_ROIcount = 100*((L_ROIcount-R_ROIcount)/(L_ROIcount+R_ROIcount));
        else
            LI_ROIcount = inf;
        end
        Thrshd_LI_ROIcount(k,2)=LI_ROIcount;
    end
    switch time_interval
        case {2; 1}
            subplot(4,4,plot_ind)
            plot_ind=plot_ind+1;
            plot(Thrshd_LI_ROIcount(:,1),Thrshd_LI_ROIcount(:,2));
            title([RoiLabels{ii} ' Count-based']);
    end
    
    % Run a loop to plot LIs based on space-time voxel-magnitude average as a function of threshold. But this metrix does not seem
    % to be very useful because the LIs tend to be around 0.
    k=0;
    %JL@10/30/14   Rng= threshold:1:AllMax;
    for thrTmp = Rng
        ind = LHvals > thrTmp;
        LHvals_aboveThreshold = LHvals(ind); % a 1-D matrix, no need for mean(mean()) later
        ind = RHvals > thrTmp;
        RHvals_aboveThreshold = RHvals(ind); % a 1-D matrix, no need for mean(mean()) later
        if isempty(LHvals_aboveThreshold)
            L_ROIavg=0;  %to prevent error of "Warning: Divide by zero" when LHvals_aboveThreshold is an empty matrix and you are doing operation of mean(LHvals_aboveThreshold)
        else
            L_ROIavg=mean(LHvals_aboveThreshold);
        end
        if isempty(RHvals_aboveThreshold)
            R_ROIavg=0;  %to prevent error of "Warning: Divide by zero" when RHvals_aboveThreshold is an empty matrix and you are doing operation of mean(RHvals_aboveThreshold)
        else
            R_ROIavg=mean(RHvals_aboveThreshold);
        end
        
        k=k+1;
        Thrshd_LI_ROIavg(k,1)=thrTmp;
        if L_ROIavg+R_ROIavg ~= 0
            LI_ROIavg = 100*((L_ROIavg-R_ROIavg)/(L_ROIavg+R_ROIavg));
        else
            LI_ROIavg = inf;
        end
        Thrshd_LI_ROIavg(k,2)=LI_ROIavg;
    end
    
    switch time_interval
        case {2; 1}
            subplot(4,4,plot_ind)
            plot_ind=plot_ind+1;
            plot(Thrshd_LI_ROIavg(:,1),Thrshd_LI_ROIavg(:,2)); %JS 092815 changed plot to subplot
            title([RoiLabels{ii} ' Average-based']);
    end
    
end
set(gcf, 'Position', [500   500   1000   800]);

%%
[name, ~, ~] = fileparts(sResultP.SurfaceFile);

% Split the name into first name and last name
name_parts = split(name, '_');

% Convert the first letter of each name to uppercase
first_name = [upper(name_parts{1}(1)) , name_parts{1}(2:end)];
last_name = upper(name_parts{2}(1));

% Combine the names and format as per requirement
formatted_name = [first_name '_' last_name];

folderPath = fullfile('/MEG_data/LanguageLI',formatted_name);

disp(folderPath)
if ~exist(folderPath, 'dir')
    mkdir(folderPath);
    disp('Folder created successfully.');
else
    disp('Folder already exists.');
end

%% Write out data to excel file  JS 09/28/15
disp('enter saving file, e.g., DFNM_LCVM or PN_dSPM')
sname = input('','s');
if time_interval == 3
    tempfile=fopen(strcat(fullfile(folderPath, ['/LI_ROItable_',sname,'_thresh',num2str(Ratio4Threshold),'.xls'])),'w') ;
else
    tempfile=fopen(strcat(fullfile(folderPath,['/LI_ROItable_',sname,'_','time',num2str(timerange(1)),'-',num2str(timerange(2)),'_thresh',num2str(Ratio4Threshold),'.xls'])),'w') ;
end

fprintf(tempfile,'%s\t',LI_label_out{:});
fprintf(tempfile,'\n');
fprintf(tempfile,'%f\t',Summ_LI);
fprintf(tempfile,'\n\nThreshold');
fprintf(tempfile,'%f\t',threshold);
fclose(tempfile);
threshold

cd(folderPath)

%% Added by VZ, display LI values
disp('=================')
disp('                 ')
a = table(RoiLabels'); a.Properties.VariableNames{'Var1'} = 'ROI';
b = table(Summ_LI'); b.Properties.VariableNames{'Var1'} = 'LI';
c = table([L_count;R_count]'); c.Properties.VariableNames{'Var1'} = 'Left_vs_right';
d = [a,b,c];
disp(d)

disp('============')
disp('To edit the LI script, first ensure Brainstorm is running. Then, open process_computeLI.m in Matlab.')
disp('Pipeline update: 08/01/23')

end
