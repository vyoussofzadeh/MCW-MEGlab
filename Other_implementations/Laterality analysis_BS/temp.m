%Script for computing laterality indices from MEG language activation dSPM maps:
%The script expects that you have exported from brainstorm to matlab a full dSPM source map <cannot just be kernel; make sure you keep the orientation normal to cortical surface by
%choosing the "Constrained (Normal/cortex)" option in the "Minimum Norm Options">, and called it 'lang'

%   HISTORY
%   10-30-14    ZL          Modified the original "LI_code_lab_v3.m" to make it work for dSPM counting at a given threshold
%   11-21-14    ZL          Added/modified some lines to display a summary table of LIs-by-ROIs at the end      
%   09-21-14    MR          Copied over edits from LI_code_lab_v4.m" to use the gui to get lang file
%   06-05-19    CU          Changed Timerange for DFNM to from [300:1.19] to [550:1.19]. 

% cd /MEG_data/LanguageLI/
clc, clear
[filename, pathname]=uigetfile('Find the dSPM file'); 
[folderpath, name, ext]=fileparts(filename);
cd(pathname);

%%
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


disp('==========')
disp('1: default timing')
disp('2: define new timing')
disp('3: averaged sources')
time_interval = input('Enter time interval: ');
switch time_interval
    case 1
        disp(['[ ', num2str(timerange(1)),',', num2str(timerange(2)), '] sec was selected for LI analysis'])
    case 2
        timerange=input('Enter Timerange: eg. [.100 .900]');
    case 3
        timerange=1;  % Make sure  these are less than 1 for mS
    otherwise
        error('There is not timerange established - pick a datatype')
end
        
% if datatype==1
% elseif datatype==2
%     timerange=[.550 1.19] %%CUstine - changed on June 5th, 2019 from [300:1.19] to [550:1.19]. 
% elseif datatype==3
%     timerange=[.200 .600] %%default 200 - 600ms!
% elseif datatype==4
%     timerange=input('Enter Timerange: eg. [.100 .900]')
% else 
%     error('There is not timerange established - pick a datatype')
% end

lang=load(filename);
font_size_graphs=8;
disp('==========')

%% Check to make sure this is a dSPM file (assuming you already loaded a BST-exported structure called "lang"):
%  if(strcmp(lang.Function,'dspm')),
%      disp('Good, this is not a zscore file!'); 
%  else
%      disp('I am not expecting a zscore file but this seems to be one! Please check!'); return;
%  end
%  
%% define calculation parameters:  sampling rate, pretrigger baseline, time window of interest, etc.
switch time_interval
    case 3
%         samplerate=round(inv((lang.time(end)-lang.time(1))/length(lang.time)))-1
        samplerate = 1;
    otherwise
        %samplerate=input('Enter sampling rate of the recording in Hz:\n');
        %samplerate=round(inv((lang.Time(end)-lang.Time(1))/length(lang.Time)),-2) %JS 11/11/15
        samplerate=round(inv((lang.Time(end)-lang.Time(1))/length(lang.Time)))-1
        %pretrig=input('Enter pretrigger baseline imported into brainstorm in ms (a positive value, e.g. 200): \n');
        %startwindow=input('Enter t1 (start of time window of interest) relative to stim onset in ms: \n');
        %endwindow=input('Enter t2 (end of time window of interest) relative to stim onset in ms: \n');
end
%% Estimate the p-value corresponding to the threshold applied to Z-scores--for future estimate this from a function fit
%if threshold==5
%    pVal=2.87e-7;
%elseif threshold==4
%    pVal=0.000032;
%elseif threshold ==3.719
%    pVal=0.0001;
%elseif threshold ==3.09
%    pVal=0.001;
%else
%    pVal=NaN;
%end; 
%fprintf('\n\nLaterality indices are based on a p value of %g.\n\n', pVal);

TotROI=8;  %define the number of ROIs to be used for LI-calculation. % added JL 11212014
Ratio4Threshold= 0.5; % ratio that will be multiplied by the maximum dSPM value across all ROIs within the time window of [startwindow,endwindow], in order to get a threshold for cutting off non-significant voxels. Added JL@10/30/14.
%JL@10/30/14    p_threshold = 0.001 % 1e-10 %  0.05  % ; 1e-7 % % the threshold (with Bonferroni correction) --- voxels with p less than this (after Bonferroni correction) is significant and used for LI calculation. Added JL 07172013


%% load the dSPM data and use the absolute value for these dSPM values
ImageGridAmp = lang.ImageGridAmp;
ImageGridAmp=abs(ImageGridAmp);  % added JL 10/30/14

%% define the window of interest based on sample rate ans pretrigger data imported into brainstorm; get some max values over various time windows
switch time_interval
    case 3
        GlobalMax=max(max(ImageGridAmp)); % max value over all time points. Added JL @ 10/30/14
        AllMax = max(max(ImageGridAmp)); % max value over the time window of interest
    otherwise
        %t1 = floor((startwindow+pretrig)*samplerate/1000);%start window of interest
        t1=find(lang.Time >= timerange(1),1) % JS 11/11/15
        %t2 = floor((endwindow+pretrig)*samplerate/1000);%end window of interest
        t2=find(lang.Time >= timerange(2),1) % JS 11/11/15
        AllMax = max(max(ImageGridAmp(:,t1:t2))); % max value over the time window of interest
        GlobalMax=max(max(ImageGridAmp)); % max value over all time points. Added JL @ 10/30/14
end

%% Convert DesikenKilliany scout to select scouts
scoutFile='./scout_Desikan-Killiany_68.mat';
sScout = load(scoutFile);

%% Define Rois  Added 10/1/15 JStout
AngSmg=[15,16,63,64]
Front=[3,4,5,6,11,12,25,26,29,30,33,34,37,38,39,40,41,42,49,50,53,54,55,56,57,58]
LatFront=[5,6,11,12,37,38,39,40,41,42,55,56,57,58]
LatTemp=[1,2,17,18,31,32,61,62,65,66,67]
% LatTemp=[1,2,17,18,31,32,61,62,65,66,67,68]
PeriSyl=[15,16,37,38,41,42,61,62,63,64]
Tanaka=[37,38,41,42,61,62,63,64]
% Temp=[1,2,9,10,13,14,17,18,19,20,27,28,31,32,35,36,61,62,65,66,67,68]
Temp=[1,2,9,10,13,14,17,18,19,20,27,28,31,32,35,36,61,62,65,66,67]
% Whole=[1:68]
Whole=[1:67]

RoiLabels={'AngSmg', 'Front','LatFront','LatTemp', 'PeriSyl', 'Tanaka','Temp','Whole'}
RoiIndices={AngSmg, Front,LatFront,LatTemp, PeriSyl, Tanaka,Temp,Whole}

% RoiLabels={'AngSmg', 'Front','LatFront', 'PeriSyl', 'Tanaka'}
% RoiIndices={AngSmg, Front,LatFront, PeriSyl, Tanaka}
%%%%%%%%%%%%% To verify the accuracy of the above regions, uncomment the
%%%%%%%%%%%%% below code and it will print out the region labels to an
%%%%%%%%%%%%% excel file
% tempfile=fopen(strcat('./ROI_subtable.xls'),'w')
% %print each subregion for verification of labels
% RoiLabels={'AngSmg', 'Front','LatFront','LatTemp', 'PeriSyl', 'Tanaka','Temp','Whole'}
% RoiIndices={AngSmg, Front,LatFront,LatTemp, PeriSyl, Tanaka,Temp,Whole}
% 
% for i=1:length(RoiLabels)
%     fprintf(tempfile,'%s\n',RoiLabels{i})
%     fprintf(tempfile,'%s\t',sScout.Scouts(RoiIndices{i}).Label)
%     fprintf(tempfile,'\n')
% end
% fclose(tempfile)
%% Load the scout file(s)& compute LI 
s1='LI_';
Summ_LI=zeros(1,TotROI); % initialize the vector that summarizes the final LIs  % added JL 11212014
Summ_LI_Label='ROI Labels: '; % initialize the string that summarizes the ROI labels  % added JL 11212014
LI_label_out={}
figure
plot_ind=1;
% threshold 

for ii = 1:8 % 1:8
        s2=RoiLabels{ii};
        %Odd indices are left Rois
        Ltemp_region=[];
        Ltemp_label=[];
        hemi_roi_num=length(RoiIndices{ii})
        curr_subregion=sScout.Scouts(RoiIndices{ii})
        for i=1:2:hemi_roi_num
            Ltemp_region=[Ltemp_region,curr_subregion(i).Vertices];
            Ltemp_label=[Ltemp_label, curr_subregion(i).Label];
        end
        %Even indices are right Rois
        Rtemp_region=[];
        Rtemp_label=[]
        for i=2:2:hemi_roi_num
            Rtemp_region=[Rtemp_region,curr_subregion(i).Vertices];
            Rtemp_label=[Rtemp_label, curr_subregion(i).Label];
        end                   

        LHscout = Ltemp_region;
        RHscout = Rtemp_region;
        
        %LHscout = sScout.Scouts(1).Vertices;
        %RHscout = sScout.Scouts(2).Vertices;

        % Compute laterality index based on a probability threshold
        % choose threshold based on a criterion p value of the z-score 
        % 
        %---threshold---%
        % Threshold of 3.719 corresponds to a CDF value of ~0.9999 (i.e
        % p<0.0001), 4 corresponds to a CDF value of ~0.999968 (p<0.000032) & 3.09 
        % corresponds to a CDF value of ~0.999 (p<0.001)

        % NOTE: The positive one-sided threshold ensures we are looking at
        % areas where language activation is greater than that for the control
        % condition

        switch time_interval
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
        

        % For maps from each side find indices of above-threshold voxels
        % With positive thresholds we are rectifying the maps so only areas with language >
        % noise by threshold amount is considered in further LI calculations 
        disp('===================>')  % to separate output from different ROIs
        
        
        %JL@10/30/14 p_Bonferroni=p_threshold/sum(length(LHscout)+length(RHscout)); % Bonferroni-corrected p-threshold. Added JL 07172013
        %JL@10/30/14 threshold = norminv(1-p_Bonferroni,0,1);  % Z-threshold after Bonferroni correction to the desired p-value for this ROI. Added JL 07172013 
        
        %threshold = Ratio4Threshold*GlobalMax;  % dSPM threshold to get rid of non-significant voxels. Added JL@10/30/14
        threshold = Ratio4Threshold*AllMax;  % dSPM threshold to get rid of non-significant voxels. Added JL@10/30/14
%        threshold = Ratio4Threshold*ROIMax;
        %JL@10/30/14   fprintf('%s:  Laterality indices are based on a p value of %g, which correspond to a Z-score threshold of %g, \n after Bonferroni correction using total vertices of %g in this ROI (L: %g, R: %g).\n\n', s2, p_threshold,threshold,sum(length(LHscout)+length(RHscout)),length(LHscout),length(RHscout));
        %JL@10/30/14   fprintf('The Bonferroni-corrected p threshold becomes %g.\n\n', p_Bonferroni);
        fprintf('%s:  Laterality indices are based on a dSPM threshold of %g  using total vertices of %g in this ROI (L: %g, R: %g).\n\n', s2, threshold,sum(length(LHscout)+length(RHscout)),length(LHscout),length(RHscout)); % added JL@10/30/14

        ind_L = find(LHvals > threshold);
        ind_R = find(RHvals > threshold);  

        % ROI count --- above threshold voxels only
        %JL@10/30/14    indHalfROImaxL = find(LHvals > ROIMax/2);%indHalfROImaxL = find(LHvals > ROIMax/2); % Not being used right now
        L_ROIcount = length(ind_L);
        %JL@10/30/14    indHalfROImaxR = find(RHvals > ROIMax/2); % Not being used right now
        R_ROIcount = length(ind_R);
        ROIcount=sum(L_ROIcount+R_ROIcount); % to report total significant voxels over space-time 
        LI_ROIcount = 100*((L_ROIcount-R_ROIcount)/(L_ROIcount+R_ROIcount));
        eval([strcat(s1,s2) '=' 'LI_ROIcount'])
        fprintf('   Based on count of %d significant voxels over space & time on the R and %d on the L\n',R_ROIcount,L_ROIcount);
        Summ_LI(ii)=LI_ROIcount;  % added JL 11212014
        Summ_LI_Label=[Summ_LI_Label  sprintf('\t') s2];   % added JL 11212014
        LI_label_out=[LI_label_out, s2]
        

        % ROI average --- above threshold voxels only
        LHvals_aboveThreshold = LHvals(ind_L); % a 1-D matrix, no need for mean(mean()) later
        RHvals_aboveThreshold = RHvals(ind_R); % a 1-D matrix, no need for mean(mean()) later
        L_ROIavg=mean(LHvals_aboveThreshold);
        R_ROIavg=mean(RHvals_aboveThreshold);
        ROIavg = mean([L_ROIavg,R_ROIavg]);
        LI_ROIavg = 100*((L_ROIavg-R_ROIavg)/(L_ROIavg+R_ROIavg));
        eval([strcat(s1,s2) '=' 'LI_ROIavg'])
        fprintf('   Based on average of %d significant voxels over space & time on the R and %d on the L\n',R_ROIcount,L_ROIcount);


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
        %figure;
        subplot(4,4,plot_ind)
        plot_ind=plot_ind+1;
        plot(Thrshd_LI_ROIcount(:,1),Thrshd_LI_ROIcount(:,2));
        %xlabel(['Threshold of dSPM for ROI ' s2 '. ROIMax is ' num2str(ROIMax) ' AllMax is ' num2str(AllMax) ],'FontSize',font_size_graphs);  %xlabel(['Threshold of Z-score for ROI ' s2 '. ROIMax is ' num2str(ROIMax) ' AllMax is ' num2str(AllMax) ]); % modified JL@10/30/14
        xlabel(['ROI_Threshold:' s2 '-ROIMax is:' num2str(ROIMax) '-AllMax is:' num2str(AllMax) ],'FontSize',font_size_graphs);  %xlabel(['Threshold of Z-score for ROI ' s2 '. ROIMax is ' num2str(ROIMax) ' AllMax is ' num2str(AllMax) ]); % modified JL@10/30/14
        %ylabel(['Count-based LI using ' scoutFile],'FontSize',font_size_graphs);
        ylabel([scoutFile],'FontSize',font_size_graphs);
        title([RoiLabels{ii} ' Count-based']);

        % Run a loop to plot LIs based on space-time voxel-magnitude average as a function of threshold. But this metrix does not seem
        % to be very useful because the LIs tend to be around 0.
        k=0;
        %JL@10/30/14   Rng= threshold:1:AllMax; 
        for thrTmp = Rng,
            ind = find(LHvals > thrTmp);
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
        %figure;
        subplot(4,4,plot_ind)
        plot_ind=plot_ind+1;
        plot(Thrshd_LI_ROIavg(:,1),Thrshd_LI_ROIavg(:,2)); %JS 092815 changed plot to subplot
        xlabel(['Threshold of Z-score for ROI ' s2 '. ROIMax is ' num2str(ROIMax) ' AllMax is ' num2str(AllMax) ],'FontSize',font_size_graphs);
        ylabel(['Average-based LI using ' scoutFile],'FontSize',font_size_graphs);
        title([RoiLabels{ii} ' Average-based'])
        
        
end

disp([sprintf('\n\n') Summ_LI_Label]);   % added JL 11212014
Summ_LI_text='Count-based LIs:' ;  %initialize the text string that summarizes the LIs for ROIs % added JL 11212014
Summ_LI_textParen='Count-based LIs:' ;  %initialize the text string that summarizes the LIs for ROIs w/ parenthesis % added JL 11212014
for jj =  1:TotROI  % added JL 11212014
    Summ_LI_text=[Summ_LI_text  sprintf('\t') num2str(Summ_LI(jj),'%10.1f')];  
end
disp(Summ_LI_text)  % added JL 11212014

%for kk =  1:TotROI  % added JL 11212014
%    Summ_LI_textParen=[Summ_LI_textParen  sprintf('\t(') num2str(Summ_LI(kk),'%10.1f') ')'];  
%end
disp(Summ_LI_textParen)  % added JL 11212014

%% Write out data to excel file  JS 09/28/15
tempfile=fopen(strcat('./ROI_table_',name,'.xls'),'w') ;
fprintf(tempfile,'%s\t',LI_label_out{:});
fprintf(tempfile,'\n');
fprintf(tempfile,'%f\t',Summ_LI);
fprintf(tempfile,'\n\nThreshold');
fprintf(tempfile,'%f\t',threshold);
fclose(tempfile);
threshold



