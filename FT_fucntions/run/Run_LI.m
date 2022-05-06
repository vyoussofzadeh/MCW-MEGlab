
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
effect = 1;

% lang = load(filename);
font_size_graphs=8;
% disp('==========')

samplerate = 1;

TotROI=8;  %define the number of ROIs to be used for LI-calculation. % added JL 11212014
Ratio4Threshold= 0.60; % ratio that will be multiplied by the maximum dSPM value across all ROIs within the time window of [startwindow,endwindow], in order to get a threshold for cutting off non-significant voxels. Added JL@10/30/14.

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
AngSmg   = [15,16,63,64];
Front    = [3,4,5,6,11,12,25,26,29,30,33,34,37,38,39,40,41,42,49,50,53,54,55,56,57,58];
LatFront = [5,6,11,12,37,38,39,40,41,42,55,56,57,58];
LatTemp  = [1,2,17,18,31,32,61,62,65,66,67,68];
PeriSyl  = [15,16,37,38,41,42,61,62,63,64];
Tanaka   = [37,38,41,42,61,62,63,64];
Temp     = [1,2,9,10,13,14,17,18,19,20,27,28,31,32,35,36,61,62,65,66,67,68];
Whole    = 1:68;
Selective = [37,38,45,46,11,12,65,66];
% Selective2 = [17,18,31,32,61,62,65,66];

% Selective = [33,34,37,38,41,42,61,62,63,64];
% Selective = [35,36,37,38,41,42,61,62,63,64];
% Selective = [13, 14, 49, 50, 47,48, 35, 36, 67, 68];

RoiLabels = {'AngSmg', 'Front','LatFront','LatTemp', 'PeriSyl', 'Tanaka','Temp','Whole','Selective'};
RoiIndices = {AngSmg, Front,LatFront,LatTemp, PeriSyl, Tanaka, Temp, Whole, Selective};

%% Load the scout file(s)& compute LI
s1='LI_';
Summ_LI=zeros(1,TotROI); % initialize the vector that summarizes the final LIs  % added JL 11212014
Summ_LI_Label='ROI Labels: '; % initialize the string that summarizes the ROI labels  % added JL 11212014
LI_label_out={};
% figure
plot_ind=1;
% threshold

for ii = 1:length(RoiLabels)
    
    s2 = RoiLabels{ii};
    %Odd indices are left Rois
    Ltemp_region = [];
    Ltemp_label  = [];
    hemi_roi_num=length(RoiIndices{ii});
    curr_subregion=sScout.Scouts(RoiIndices{ii});
    
    %--- Added by vyz
    region_all = []; region_mean = []; label = [];
    for i=1:hemi_roi_num
        region_all{i} = ImageGridAmp(curr_subregion(i).Vertices);
        region_mean(i) = mean(ImageGridAmp(curr_subregion(i).Vertices));
        label{i} = curr_subregion(i).Label;
    end
    
    parcel = [];
    parcel.label = label;
    parcel.val = region_all;
    parcel.valmean = region_mean;
    parcel.atlas = RoiLabels{ii};
    parcel_all{ii} = parcel;
    
    %---
    k = 1;
    for i=1:2:hemi_roi_num
        Ltemp_region = [Ltemp_region,curr_subregion(i).Vertices];
        Ltemp_label{k} = curr_subregion(i).Label; k = k+1;
    end
    
    %Even indices are right Rois
    Rtemp_region = [];
    Rtemp_label  = [];
    
    k = 1;
    for i=2:2:hemi_roi_num
        Rtemp_region = [Rtemp_region,curr_subregion(i).Vertices];
        Rtemp_label{k} = curr_subregion(i).Label; k = k+1;
    end
    LHscout = Ltemp_region;
    RHscout = Rtemp_region;
    
    %First parse the maps into separate space-times maps for each side
    LHvals = ImageGridAmp(LHscout);
    LH_max = max(max(LHvals));
    RHvals = ImageGridAmp(RHscout);
    RH_max = max(max(RHvals));
    ROIMax = max(LH_max,RH_max);
    
    
    threshold = Ratio4Threshold*AllMax;  % dSPM threshold to get rid of non-significant voxels. Added JL@10/30/14
    
    ind_L = find(LHvals > threshold);
    ind_R = find(RHvals > threshold);
    
    L_ROIcount = length(ind_L); L_count(ii) = L_ROIcount;
    R_ROIcount = length(ind_R); R_count(ii) = R_ROIcount;
    ROIcount=sum(L_ROIcount+R_ROIcount); % to report total significant voxels over space-time
    LI_ROIcount = 100*((L_ROIcount-R_ROIcount)/(L_ROIcount+R_ROIcount));
    Summ_LI(ii) = LI_ROIcount;  % added JL 11212014
    Summ_LI_Label = [Summ_LI_Label  sprintf('\t') s2];   % added JL 11212014
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
    
    k=0;
    %JL@10/30/14   Rng= threshold:1:AllMax;
    for thrTmp = Rng
        ind = LHvals > thrTmp;
        LHvals_aboveThreshold = LHvals(ind); % a 1-D matrix, no need for mean(mean()) later
        ind = find(RHvals > thrTmp);
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
    
end

%% Write out data to excel file  JS 09/28/15
% tempfile=fopen(strcat('./ROI_table_',name,'.xls'),'w') ;
% fprintf(tempfile,'%s\t',LI_label_out{:});
% fprintf(tempfile,'\n');
% fprintf(tempfile,'%f\t',Summ_LI);
% fprintf(tempfile,'\n\nThreshold');
% fprintf(tempfile,'%f\t',threshold);
% fclose(tempfile);

%%
disp('                 ')
a = table(RoiLabels'); a.Properties.VariableNames{'Var1'} = 'ROI';
b = table(Summ_LI'); b.Properties.VariableNames{'Var1'} = 'LI';
c = table([L_count;R_count]'); c.Properties.VariableNames{'Var1'} = 'Left_vs_right';
d = [a,b,c];
% disp(d)

%%
% for i=1:length(parcel_all)
%     figure, bar(parcel_all{i}.valmean);
%     L = length(parcel_all{i}.valmean);
%     set(gca,'Xtick', 1:L,'XtickLabel',parcel_all{i}.label);
%     box off
%     set(gca,'color','none');
%     xlim([0,L+1])
%     xlabel('ROI');
%     ylabel('Source power');
%     set(gcf, 'Position', [100   100   1500   500]);
%     set(gca,'FontSize',10,'XTickLabelRotation',90);
%     grid
%     title(parcel_all{i}.atlas)
% end
