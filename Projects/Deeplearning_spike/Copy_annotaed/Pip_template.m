%% The Spike Detection MEG pipline

% Spike detection MEG pipline
% Writtern by MCW group, Youssofzadeh, Vahab <vyoussofzadeh@mcw.edu>
% Lastest update: 08/09/2022

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
flag.analysis = 1;

%% Datalog (subject details)
Datalog = [];

%% Initial settings
cd '/MEG_data/LAB_MEMBERS/Vahab/Github/MCW-MEGlab/FT';
restoredefaultpath
cd_org = cd;
addpath(genpath(cd_org));

%- Input dir
indir = '/MEG_data/Research_studies/Epil_annotated_data/annotated_info';
%- Output dir
outdir = '/MEG_data/Research_studies/Epil_annotated_data/annotated_data_nospike';

%- Adding path
cfg_init = [];
cfg_init.path_tools = '/MEG_data/LAB_MEMBERS/Vahab/Github/tools';
[allpath, atlas] = vy_init(cfg_init);

addpath('/MEG_data/LAB_MEMBERS/Vahab/Github/MCW-MEGlab/MCW_MEGlab_git/external')

%%
datadir = '/MEG_data/Research_studies/Epil_annotated_data/annotated_data';
cd(datadir)
d = rdir([datadir,'/*.mat']);

%%
clear subj run sub_run
for i=1:length(d)
    [pathstr, name] = fileparts(d(i).name);
    tkz = tokenize(name,'_');
    subj{i} = [tkz{1}, '_', tkz{2}];
    sub_run{i,:} = [subj{i}, '_', tkz{3}];
end
[sub_run_unq,IA,IC] = unique(sub_run);
% disp(sub_run_unq);

sub_run_unq1 = [];
for i=1:length(sub_run_unq)
    sub_run_unq1{i} = [num2str(i), '_', sub_run_unq{i}];
end
disp(sub_run_unq1')

%%
% disp('Enter sub')
% subid = input('');

%%
avg_eeg = zeros(1,68);
avg_meg = zeros(1,68);

pca_eeg_sub = [];
pca_meg_sub = [];

for i= 1:length(d)
    disp([num2str(i),'/',num2str(length(d))])
    [pathstr, name] = fileparts(d(i).name);
    load(d(i).name);
    
    for j=1:length(anot_data_all)
        anot_data = anot_data_all{j};
        
        cfg = [];
        cfg.channel = 'EEG*';
        eeg = ft_selectdata(cfg, anot_data);
        %         figure, plot(smooth(mean(eeg.trial{:},1)))
        
        cfg = [];
        cfg.channel = 'MEG*';
        meg = ft_selectdata(cfg, anot_data);
        %         figure, plot(smooth(mean(meg.trial{:},1)))
        
        avg_eeg(j,:) = zeros(1,68); D_eeg = smooth(mean(eeg.trial{:},1)); avg_eeg(j,1:length(D_eeg)) = D_eeg;
        avg_meg(j,:) = zeros(1,68); D_meg = smooth(mean(meg.trial{:},1)); avg_meg(j,1:length(D_meg)) = D_meg;
        
%         pca_eeg(j,:) = smooth(smooth(do_pca(eeg.trial{:}, 1)));
%         pca_meg(j,:) = smooth(smooth(do_pca(meg.trial{:}, 1)));

        
        %         cfg            = [];
        %         cfg.method     = 'runica';
        %         cfg.numcomponent = 5;       % specify the component(s) that should be plotted
        %         comp           = ft_componentanalysis(cfg, eeg);
        %         figure, plot(smooth(abs(comp.trial{1}(1,:))))
        
        %         cfg = [];
        %         cfg.blocksize = anot_data.time{1}(end) - anot_data.time{1}(1);
        %         cfg.viewmode = 'vertical'; %butterfly';
        %         cfg.continuous = 'yes';
        %         cfg.axisfontsize = 7;
        %         cfg.fontsize = 7;
        %         cfg.channel = 'EEG*';
        %         cfg.preproc.demean = 'yes';
        %         cfg.position = [300   900   500   1500];
        %         ft_databrowser(cfg, anot_data);
        %         cfg.channel = 'MEG*';
        %         cfg.position = [850   900   500   1500];
        %         ft_databrowser(cfg, anot_data);
        
        %         pause,
%         close all,
    end
%     avg_eeg_sub {i} = avg_eeg(j,:);
%     avg_meg_sub {i} = avg_meg(j,:);
    
    pca_eeg_sub {i} = pca_eeg(j,:);
    pca_meg_sub {i} = pca_meg(j,:);
end

% figure,plot(mean(abs(avg_eeg(1:67)),1))
% figure,plot(mean(abs(avg_meg(1:67)),1))

%%
sdir = '/MEG_data/Research_studies/Epil_annotated_data/datasave';
save(fullfile(sdir,'annotated_data'),'avg_eeg_sub','avg_meg_sub')
% load((fullfile(sdir,'annotated_data')))

%%
for i=1:length(avg_eeg_sub)

    tmp = avg_eeg_sub{1};
        pca_val = smooth(smooth(do_pca(eeg.trial{:}, 1)));
        figure, plot(pca_val)
end

%%
addpath('/MEG_data/LAB_MEMBERS/Vahab/Github/tools/megclinic_development/tools/func')

% for i=1:length(avg_eeg_sub)
%     pca_val = do_pca(avg_eeg_sub{i}, 1);
% end
a = load('sample_pca.mat');
A = interp1(1:length(pca_val),pca_val,1:length(a.comp));
figure, plot(y2i)
figure,plot(mean(a.comp,1))

B = mean(a.comp,1);
indx = ~(isnan(A) | isnan(B));
r = corr2(abs(A(indx)),abs(B(indx)));

%%
numcomponent = 1;
Nchans = size(avg_eeg,2);
dat = avg_eeg';
C = (dat*dat')./(size(dat,2)-1);

% eigenvalue decomposition (EVD)
[E,D] = eig(C);

% sort eigenvectors in descending order of eigenvalues
d = cat(2,(1:1:Nchans)',diag(D));
d = sortrows(d, -2);

% return the desired number of principal components
unmixing = E(:,d(1:numcomponent,1))';
mixing = [];

for j=1:numcomponent
    figure, plot(smooth(unmixing(j,1:67)))
end


%%
numcomponent = 1;
Nchans = size(avg_meg,2);
dat = avg_meg';
C = (dat*dat')./(size(dat,2)-1);

% eigenvalue decomposition (EVD)
[E,D] = eig(C);

% sort eigenvectors in descending order of eigenvalues
d = cat(2,(1:1:Nchans)',diag(D));
d = sortrows(d, -2);

% return the desired number of principal components
unmixing = E(:,d(1:numcomponent,1))';
mixing = [];

for j=1:numcomponent
    figure, plot(smooth(unmixing(j,1:67)))
end

%%

% [U,S,V] = svd(avg_eeg);
% [U,S,V] = svd(avg_eeg,"econ");
%
% comp = pca(avg_eeg);
%
% figure, plot(s)


%%