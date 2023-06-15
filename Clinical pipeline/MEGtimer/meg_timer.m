% Timing info of the CNRP MEG data
% Writtern by MCW group, Youssofzadeh, Vahab <vyoussofzadeh@mcw.edu>
% Lastest update: 12/13/2022

%%
addpath('/MEG_data/LAB_MEMBERS/Vahab/Github/MCW-MEGlab/MCW_MEGlab_git/FT_fucntions/functions/External')
addpath('/usr/local/MATLAB_Tools/mne')

cd_datapath = '/MEG_acq/cnrp_tacs_healthy';
% cd_datapath = '/MEG_acq/cnrp_tacs_healthy/hc004_v3';
% cd_datapath = '/MEG_acq/cnrp_tacs_healthy/hc008_v1';
% 
% cd_datapath = '/MEG_acq/cnrp_tacs_healthy/mcwa023_v1';
% cd_datapath = '/MEG_acq/cnrp_tacs_healthy/mcwa023_v2';

%% HC data
cd(cd_datapath)
dd = dir('./hc*');
for j = 1:length(dd)
    cd(cd_datapath)
    cd(dd(j).name)
    d = rdir('./**/tsss/*.fif');
    tt = []; date_data = [];
    for i=1:length(d)
        data = fiff_setup_read_raw(d(i).name);
        tt(i) = double(data.last_samp)/data.info.sfreq;
        date_data{i} = d(i).date(1:12);
    end    
    summ_tim_hc (j) = sum(tt)/60 + length(d);
    data_name_hc{j} = dd(j).name;
    date_data_hc{j} = unique(date_data);
end

%% MCW (pataint) data
cd(cd_datapath)
dd = dir('./mcw*');
for j = 1:length(dd)
    cd(cd_datapath)
    cd(dd(j).name)
    d = rdir('./**/tsss/*.fif');
    tt = [];  date_data = [];
    for i=1:length(d)
        data = fiff_setup_read_raw(d(i).name);
        tt(i) = double(data.last_samp)/data.info.sfreq;
        date_data{i} = d(i).date(1:12);
        disp(d(i).date)
    end    
    summ_tim_mcw (j) = sum(tt)/60 + length(d);
    data_name_mcw{j} = dd(j).name;
    tmp = unique(date_data);
    date_data_mcw{j} = tmp(1);
end

%% Merging details
clc
summ_tim  = [summ_tim_mcw, summ_tim_hc];
data_name = [data_name_mcw, data_name_hc];
date_all = [date_data_mcw, date_data_hc];

timehr = (summ_tim/60)';

for i = 1:length(timehr)
    if timehr(i) > 1.15
        rtimehr(i) = 2;
    elseif  timehr(i) > 0.15 &&  timehr(i) < 1
        rtimehr(i) = 1;
    else
        rtimehr(i) = round(timehr(i));
    end
end

disp(data_name');
disp(summ_tim');
disp(timehr)
disp(rtimehr);
% disp(date_all)
