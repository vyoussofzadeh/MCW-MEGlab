clear; clc, close('all'); warning off

%% Flags
flag.preprocessing.filtering = 1;
flag.preprocessing.artifact = 1;
flag.preprocessing.ica = 1;
flag.notch = 1;
flag.freq = 1;     % TFR & FFT
flag.time = 1;     % Time-locked & Cov estimation
flag.gave = 0;     % grand average analysis
flag.anatomy = 1;     % grand average analysis
flag.sourceanalysis = 1;     % grand average analysis
flag.speechanalysis = 1;     % speech analysis

%% Datalog (subject details)
Datalog = [];

%% Initial settings
% set(0,'DefaultFigureWindowStyle','docked')
% set(0,'DefaultFigureWindowStyle','normal')

cd '/MEG_data/Vahab/Github/MCW-MEGlab/FT';
restoredefaultpath
cd_org = cd;
addpath(genpath(cd_org));

%- Input dir
indir = '/MEG_data/epilepsy';
%- Output dir
outdir = '/MEG_data/Vahab/Processed_data';

%- Adding path
cfg_init = [];
cfg_init.path_tools = '/MEG_data/Vahab/Github/tools';
[allpath, atlas] = vy_init(cfg_init);

%%
cd(indir)
[subjdir] = uigetdir;

%%
d1 = rdir([subjdir,['/**/','sss','/*/*_sss.fif']]);
% d1 = rdir([subjdir,['/**/','sss','/*/*_sss*.fif']]);

for i=1:length(d1)
    alldata{i} = d1(i).name;
end
disp(alldata')

%%
disp('1: RestEC'); % Eyes-open
disp('2: RestEO'); % Eyes-closed
disp('3: Spont');
tag_sel = input('Enter data type: ');
switch tag_sel
    case 1
        tag = 'RestEC';
    case 2
        tag = 'RestEO';
    case 3
        tag = 'spont';
end

%%
d1 = rdir([subjdir,['/**/','sss','/*',tag,'*/*sss.fif']]);
d = rdir([subjdir,['/**/','sss','/*',tag,'*/*_raw.fif']]);
d2 = [d;d1];
d = d2;

%%
clear subj datafolder datafile datafile1
for i=1:length(d)
    [pathstr, name] = fileparts(d(i).name);
    datafolder{i} = pathstr;
    datafile{i} = d(i).name;
    Index = strfind(datafile{i}, '/');
    subj = datafile{i}(Index(3)+1:Index(4)-1);
end
datafile1 = datafile';

datafile2 = [];
for j=1:size(datafile1,1)
    datafile2{j} = [num2str(j), ':', '..',datafile1{j,1}(35:end)];
end
disp(datafile2')
if length(datafile1) > 1
    datasel = input('choose which data to analyze, row number:');
else
    datasel = 1;
end
disp([subj, ' and,'])
disp([datafile1{datasel}, 'was selected for the analysis ...'])
disp('============');

%%
epoch_type = 'STI101';

%% 4D layout
cfg = [];
cfg.layout = 'neuromag306mag.lay';
lay = ft_prepare_layout(cfg);
% ft_layoutplot(cfg);
disp('============');

%%
close all
datafile = datafile1{datasel}; % spm_select(inf,'dir','Select MEG folder'); % e.g. H:\VNS\MEG\C-105\CRM\1
% datafile = '/MEG_data/epilepsy/kaja_john/090701/Run_01_spontaneous_eyesclosed_raw_tsss.fif'

Index = strfind(datafile, '/');
Date  = datafile(Index(4)+1:Index(5)-1);
disp('============');
disp(datafile)
disp(['subj:',subj])
disp(['Date:',Date])
disp('============');

idx = strfind(datafile,'Run');
if ~isempty(idx)
    run  = str2double(datafile(idx(1)+3:idx(1)+4));
    run = num2str(run);
else
    idx = strfind(datafile,'run');
    if ~isempty(idx)
        run  = str2double(datafile(idx(1)+3:idx(1)+3));
        run = num2str(run);
    end
end

%%
disp(datafile)
disp(['run:',run])
disp(['subj:',subj])
disp(['Date:',Date])
disp('============');

disp('1: Yes');
disp('2: No');
inc = input('Are these correct? ');
if inc == 2
    run = input('Enter run number? ');
    subj = input('Enter subject name? ');
end

%%
%-elec/grad
sens = ft_read_sens(datafile);
sens = ft_convert_units(sens,'mm');

%%
outd.sub = fullfile(outdir,'ft_process',subj, tag, num2str(run));
if exist(outd.sub, 'file') == 0
    mkdir(outd.sub);   %create the directory
end
cd(outd.sub)
disp(['outputdir:',outd.sub])
disp('============');

%% Updating datalog
Datalog.name = subj;
Datalog.run = run;
Datalog.datafile = datafile;
Datalog.subjdir = subjdir;
Datalog.date = Date;
Datalog.outdir = outdir;

%% Preprocesssing
clear cln_data
ic_selection = 1; % 1: manual, 2: automated
Run_preprocess_SZ

%% cutting 10 sec from the end of data, button press artifacts)
disp('1: Yes');
disp('2: No');
% ct = input('Trimming data (due to button press artifacts)?');
ct = 1;
if ct==1
    %     if cln_data.time{:}(end) > 600 % 10min data
    cfg = [];
    cfg.toilim = [cln_data.time{:}(10*cln_data.fsample),cln_data.time{:}(end-10*cln_data.fsample)];
    cln_data = ft_redefinetrial(cfg,cln_data);
    %     end
end
%     ft_selectdata

%%
cfg = [];
cfg.channel = 'MEG';
cfg.covariance = 'yes';
cov_matrix = ft_timelockanalysis(cfg, cln_data);

[u,s,v] = svd(cov_matrix.cov);
%         figure;
%         semilogy(diag(s),'o-');

%%
outputmridir = fullfile(outdir,'ft_process', subj,'anat'); % output dir
if exist(outputmridir, 'file') == 0, mkdir(outputmridir); end

flag.anatomycheck = 2;
meshgridres = 2;
method = 4; % LCMV_kurtosis
Run_volumetric

if exist('./source', 'file') == 0, mkdir('./source'), end
print('./source/source_all','-dpng');

%%
array = reshape(source.avg.kurtosis, source.dim);
array(isnan(array)) = 0;
ispeak = imregionalmax(array); % findpeaksn is an alternative that does not require the image toolbox
peakindex = find(ispeak(:));
[~, i] = sort(source.avg.kurtosis(peakindex), 'descend'); % sort on the basis of kurtosis value
peakindex = peakindex(i);

npeaks = 3;
disp(source.pos(peakindex(1:npeaks),:)); % output the positions of the top peaks

pka = input('Plot peak source areas (yes=1, No=2)?');
if pka==1
    mask = 'kurtosis';
    source_int1 = source_int_all;
    tmp = abs(source_int1.(mask));
    tmp = (tmp - min(tmp(:))) ./ (max(tmp(:)) - min(tmp(:))); %
    source_int1.(mask) = tmp;
    
    for j = 1:npeaks
        cfg = [];
        cfg.funparameter = mask;
        cfg.maskparameter = mask;
        cfg.location = source.pos(peakindex(j),:); % convert from m into mm
        cfg.funcolorlim   = [0.1 1];
        cfg.opacitylim    = [0.1 1];
        %         cfg.nslices       = 16;
        %         cfg.slicerange    = [10,60];
        %         cfg.slicedim      = 3;
        cfg.opacitymap    = 'rampup';
        cfg.funcolormap =  brewermap(256, '*RdYlBu');
        ft_sourceplot(cfg, source_int1);
        print(['./source/source_all_p', num2str(j)],'-dpng');
    end
end

%% Export to AnyView format
dat = ft_fetch_data(cln_data);
hdr = ft_fetch_header(cln_data);

for i=1:size(dat,1)
    %             hdr.label{i}= ['S' num2str(i)];
    hdr.chantype{i} = 'MEG';
    hdr.chanunit{i} = 'T' ; % see note below about scaling
end

% npeaks = 5;
for i = 1:npeaks
    dat(end+1,:) = source.avg.mom{peakindex(i),:}; % see comment below about scaling
    hdr.label{end+1}= ['S' num2str(i)];
    hdr.chantype{end+1} = 'Source';
    hdr.chanunit{end+1} = 'T' ; % see note below about scaling
end
hdr.nChans = hdr.nChans+npeaks;
% ft_write_data('Case3_timeseries', dat, 'header', hdr, 'dataformat', 'anywave_ades');
dat1 = dat;

%- Marking potential spikes in the source time series
% fid = fopen('Case3_timeseries.mrk', 'w');
% fprintf(fid,'%s\r\n','// AnyWave Marker File ');
k = 1;
time_occur = [];
kk = 20;
for i = 1:1 %npeaks
    dat = source.avg.mom{peakindex(i),:};
    sd = std(dat);
    tmp = [];
    tr = zeros(size(dat));
    while length(tmp) < 5
        tr(dat>kk*sd)=1;
        [tmp, peaksample] = findpeaks(tr, 'MinPeakDistance', 300); % peaks have to be separated by 300 sample points to be treated as separate
        kk = kk-1;
        %             disp(kk);
        %             disp(peaksample)
        %             pause
    end
    for j = 1:length(peaksample)
        %         fprintf(fid, 'S%d_%02d\t', i, j); % marker name
        %         fprintf(fid, '%d\t', dat(peaksample(j))); % marker value
        %         fprintf(fid, '%d\t',source.time(peaksample(j)) ); % marker time
        %         fprintf(fid, '%d\t', 0); % marker duration
        %         fprintf(fid, 'S%d\r\n', i); % marker channel
        time_occur(j) = source.time(peaksample(j));
        k = k + 1;
        %             L = length(tmp);
        
    end
end
% fclose(fid);
disp(time_occur)

%% Sensor-level plotting
cd(outd.sub)
kk = 0.5; % Seconds
un_time_occur = unique(time_occur);

L = length(un_time_occur);
if L > 5, L =5; end

for i=1:L
    cfg = [];
    cfg.toilim = [un_time_occur(i) - kk,un_time_occur(i) + kk];
    data_spk = ft_redefinetrial(cfg, cln_data);
    
    %     kk = 0.5;
    %     cfg = [];
    %     cfg.latency = [time_occur(1) - kk,time_occur(1) + kk];
    %     data_spk2 = ft_selectdata(cfg, data_resampled);
    
    fg = [];
    cfg.blocksize = kk*2;
    cfg.viewmode = 'vertical';
    cfg.continuous = 'yes';
    ft_databrowser(cfg, data_spk);
    print(['./source/sensor_t', num2str(i)],'-dpng');
    
    %         datak = [];
    %         datak.label    = data_spk.label;
    %         datak.dimord   = 'chan';
    %         datak.kurtosis = kurtosis(data_spk.trial{1}')';
    %
    %         cfg = [];
    %         cfg.comment = 'computed channel-level kurtosis';
    %         datak = ft_annotate(cfg, datak);
    
    %         cfg = [];
    %         cfg.layout    = 'neuromag306mag.lay';
    %         cfg.parameter = 'kurtosis';
    %         figure,
    %         ft_topoplotER(cfg, datak);
    %         title([num2str(un_time_occur(i)), 'sec'])
end

%%
cd(outd.sub),
close all

%- Spikes source time series plotting
data_dummy = [];
data_dummy.trial{1} = dat1;
data_dummy.time = cln_data.time;
data_dummy.hdr = hdr;
data_dummy.label = hdr.label;

clc
tsel = [];
for j=1:L
    tsel{j} = [num2str(j), '=', num2str(un_time_occur(j))];
end
disp(tsel')
disp('1: Yes');
disp('2: No');
ETS_ask = input('Exploring a time interval and source timeseries?');


if ETS_ask ==1
    
    %     tselin = input('choose timing:');
    for j=1:length(un_time_occur)
        %
        cfg = [];
        cfg.toilim = [un_time_occur(j) - kk,un_time_occur(j) + kk];
        data_spk = ft_redefinetrial(cfg, cln_data);
        
        cfg = [];
        cfg.channel = 'MEG';
        cfg.covariance = 'yes';
        cov_matrix1 = ft_timelockanalysis(cfg, data_spk);
        
        %-
        cfg = [];
        cfg.method = 'lcmv';
        cfg.grid = individual_grid;
        cfg.headmodel = outanat.individual_headmodel;
        cfg.lcmv.keepfilter = 'yes';
        cfg.lcmv.fixedori = 'yes'; % project on axis of most variance using SVD
        cfg.lcmv.lambda = '5%';
        cfg.lcmv.kappa = 69;
        cfg.lcmv.projectmom = 'yes'; % project dipole time series in direction of maximal power (see below)
        cfg.lcmv.kurtosis = 'yes';
        source_sel = ft_sourceanalysis(cfg, cov_matrix1);
        
        cfg = [];
        cfg.mask = 'kurtosis';
        cfg.loc = 'max';
        cfg.template = outanat.mri_realigned;
        cfg.savefile = [];
        cfg.volnorm     = 2; % yes: 1
        cfg.method  = 'ortho';
        source_int_indv = vy_source_plot(cfg, source_sel);
        print(['./source/source_sel_',tsel{j}(3:end)],'-dpng');
        
        
        kk = 0.5;
        cfg = [];
        cfg.latency = [un_time_occur(j) - kk,un_time_occur(j) + kk];
        cfg.channel  = 'S*';
        data_spk1 = ft_selectdata(cfg, data_dummy);
        
        fg = [];
        cfg.blocksize = kk*2;
        cfg.channel  = 'S*';
        cfg.viewmode = 'vertical';
        cfg.continuous = 'yes';
        cfg.ylim = 'maxabs';
        ft_databrowser(cfg, data_spk1);
        title([num2str(un_time_occur(j)), 'sec']);
%         print('./source/source_timeseries','-dpng');
        print(['./source/source_timeseries_',tsel{j}(3:end)],'-dpng');
        
        %-
        datak = [];
        datak.label    = data_spk.label;
        datak.dimord   = 'chan';
        datak.kurtosis = kurtosis(data_spk.trial{1}')';
        
        cfg = [];
        cfg.comment = 'computed channel-level kurtosis';
        datak = ft_annotate(cfg, datak);
        
        cfg = [];
        cfg.layout    = 'neuromag306mag.lay';
        cfg.parameter = 'kurtosis';
        figure,
        ft_topoplotER(cfg, datak);
        title([num2str(un_time_occur(j)), 'sec'])
%         print('./source/topo_map','-dpng');
        print(['./source/topo_map_',tsel{j}(3:end)],'-dpng');
    end
end

%- Saving data,
%     array = reshape(source_sel.avg.kurtosis, source_sel.dim);
%     array(isnan(array)) = 0;
%     ispeak = imregionalmax(array); % findpeaksn is an alternative that does not require the image toolbox
%     peakindex = find(ispeak(:));
%     [peakval, i] = sort(source_sel.avg.kurtosis(peakindex), 'descend'); % sort on the basis of kurtosis value
%     peakindex = peakindex(i);
%
%     npeaks = 3;
%     disp(source_sel.pos(peakindex(1:npeaks),:)); % output the positions of the top peaks

%     pka = input('Plot peak source areas (yes=1, No=2)?');
%     if pka==1
%         mask = 'kurtosis';
%         source_int1 = source_int_indv;
%         tmp = abs(source_int1.(mask));
%         tmp = (tmp - min(tmp(:))) ./ (max(tmp(:)) - min(tmp(:))); %
%         source_int1.(mask) = tmp;
%
%         for j = 1:npeaks
%             cfg = [];
%             cfg.funparameter = mask;
%             cfg.maskparameter = mask;
%             cfg.location = source_sel.pos(peakindex(j),:); % convert from m into mm
%             cfg.funcolorlim   = [0.1 1];
%             cfg.opacitylim    = [0.1 1];
%             %         cfg.nslices       = 16;
%             %         cfg.slicerange    = [10,60];
%             %         cfg.slicedim      = 3;
%             cfg.opacitymap    = 'rampup';
%             cfg.funcolormap =  brewermap(256, '*RdYlBu');
%             ft_sourceplot(cfg, source_int1);
%             print(['./source/source_sel_p', num2str(j)],'-dpng');
%         end
%     end




%%
disp('1: Yes');
disp('2: No');
savingevent_ask = input('Saving events (readable by mnebrowse)?');
if savingevent_ask ==1
    
    %- Export events to MNEbrowse format
    cd(outd.sub)
    if exist('./event', 'file') == 0, mkdir('./event'), end
    
    hdr = ft_read_header(datafile);
    fs = hdr.orig.sfreq;
    first_samp = double(hdr.orig.raw.first_samp);
    
    clear event
    event = [first_samp, first_samp/fs, 0, 0];
    for j=1:L
        event(j+1,1) = un_time_occur(j)*fs + event(1,1);
        event(j+1,2) = un_time_occur(j) + event(1,2);
        event(j+1,3) = 0;
        event(j+1,4) = 5555;
    end
    
    tevent = table(event);
    teventcell = table2cell(tevent);
    
    savefile  = ['Evt_KurtBF_',subj,'_run', num2str(run), '_date', date,'.txt'];
    
    textfile = ['event/',savefile];
    dlmwrite(textfile, teventcell, 'delimiter','\t','precision',13);
    
    %-
    pathstr = datafolder{datasel};
    cd(pathstr)
    dlmwrite(savefile, teventcell, 'delimiter','\t','precision',13);
    
    clc
    disp('Events for data,')
    disp(datafile)
    disp('was saved at,')
    disp(pathstr)
    disp('as')
    disp(savefile)
    
    %- updating datalog,
    Datalog.event  = event;
    Datalog.time_occur = un_time_occur;
    
    %- mne-browse
    tkz = tokenize(datafile,'/');
    command = ['mbrowse ', (tkz{end})];
%     disp(tsel)
    system(command)
end

%% Saving Data-log,
cd(outd.sub)
save(fullfile(outd.sub,'Datalog.mat'), 'Datalog')

%%
%     source_temp = source_sel;
%     source_temp.pos     = template_grid.pos;
%     source_temp.dim     = template_grid.dim;
%     source_temp.inside  = template_grid.inside;
%     template_mri = ft_read_mri(fullfile(allpath.ft_path,'template/anatomy','single_subj_T1.nii')); %
%
%     cfg = [];
%     cfg.mask = 'kurtosis';
%     cfg.loc = 'max';
%     cfg.template = template_mri;
%     cfg.savefile = [];
%     cfg.volnorm     = 2; % yes: 1
%     source_int_temp = vy_source_plot(cfg, source_temp);
%
%     cfg = [];
%     cfg.subj = subj;
%     cfg.mask = 'kurtosis';
%     cfg.thre = 0.85;
%     cfg.savepath = [];
%     vy_mapvisualisation(cfg, source_int_temp);

%% DICS analysis (incomplete)
% cfg = [];
% cfg.savefile = [];
% cfg.saveflag = 2;
% cfg.foilim = [10 13];
% cfg.plotflag  = 2;
% cfg.tapsmofrq       = 1;
% cfg.taper    = 'hanning';
% f_data = vy_fft(cfg, data_spk);
% f_data.elec = sens;
%
% cfg = [];
% cfg.method = 'dics';
% cfg.grid = individual_grid;
% cfg.headmodel = outanat.individual_headmodel;
% cfg.dics.keepfilter = 'yes';
% cfg.dics.fixedori = 'yes'; % project on axis of most variance using SVD
% cfg.dics.lambda = '5%';
% cfg.dics.kappa = 69;
% % cfg.lcmv.projectmom = 'yes'; % project dipole time series in direction of maximal power (see below)
% % cfg.lcmv.kurtosis = 'yes';
% source_freq = ft_sourceanalysis(cfg, f_data);
%
%
% cfg = [];
% cfg.mask = 'pow';
% cfg.loc = 'max';
% cfg.template = outanat.mri_realigned;
% cfg.savefile = [];
% cfg.volnorm     = 2; % yes: 1
% cfg.method  = 'ortho';
% source_int_indv_freq = vy_source_plot(cfg, source_freq);

