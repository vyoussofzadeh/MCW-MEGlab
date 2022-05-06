
% cd /MEG_data/Vahab/Projects/BCI
% clear, clc, close all
% [filename, pathname]=uigetfile('Select source file');
% [folderpath, name, ext]=fileparts(filename);
% sScout = load ('./scout_Desikan-Killiany_68_template.mat');
% cd(pathname);

%%
% clc
timerange=1;  % Make sure  these are less than 1 for mS

% disp('==========')
% disp('1: Positve values');
% disp('2: Negative values');
% disp('3: Absoulte values');
% effect = input('Effects?');
effect = 3;

% lang = load(filename);
font_size_graphs=8;
% disp('==========')

samplerate = 1;

TotROI=1;  %define the number of ROIs to be used for LI-calculation. % added JL 11212014
% Ratio4Threshold= 0.60; % ratio that will be multiplied by the maximum dSPM value across all ROIs within the time window of [startwindow,endwindow], in order to get a threshold for cutting off non-significant voxels. Added JL@10/30/14.

%% load the dSPM data and use the absolute value for these dSPM values
ImageGridAmp = tmp.ImageGridAmp;

%%
switch effect
    case 1
        ImageGridAmp = ImageGridAmp;
    case 2
        ImageGridAmp = -ImageGridAmp;
    case 3
        ImageGridAmp= abs(ImageGridAmp);
end
% Edited by VY 8/7/19, absolute was commented to be account for contrating effects.
% ImageGridAmp=abs(ImageGridAmp);  % added JL 10/30/14

%% define the window of interest based on sample rate ans pretrigger data imported into brainstorm; get some max values over various time windows
GlobalMax=max(ImageGridAmp(:)); % max value over all time points. Added JL @ 10/30/14
AllMax = max(ImageGridAmp(:)); % max value over the time window of interest

%% Convert DesikenKilliany scout to select scouts
% scoutFile='./scout_Desikan-Killiany_68.mat';
% sScout = load(scoutFile);

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
for i=1:length(sScout.Scouts)
    mdk_roi(i) = mean(ImageGridAmp(sScout.Scouts(i).Vertices));
end

for i=1:l
    rois{i,:} = sScout.Scouts(i).Label;
end


