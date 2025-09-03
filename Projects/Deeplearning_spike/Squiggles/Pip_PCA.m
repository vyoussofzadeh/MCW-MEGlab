%% The Spike Detection MEG pipline

% Spike Detection MEG pipline, Deep learning process
% Writtern by MCW group, Youssofzadeh, Vahab <vyoussofzadeh@mcw.edu>
% Lastest update: 08/16/2022

clear; clc, close('all'); warning off

%% FieldTrip toolbox
restoredefaultpath % reset the default path
ft_path ='/opt/matlab_toolboxes/ft_packages/Stable_version/fieldtrip-master';
addpath(ft_path);
ft_defaults

addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/FT_fucntions/functions_new/')
addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/FT_fucntions/helper')




% restoredefaultpath % reset the default path
% ft_path ='/opt/matlab_toolboxes/ft_packages/latest/fieldtrip-master';
% addpath(ft_path);
% ft_defaults
%
% addpath('/data/MEG/Research/awang/Scripts/func')
%
datadir = '/data/MEG/Research/SpikeDectection/Epil_annotated_data/annotated_data_anonymized';

% if exist(savedir, 'file') == 0, mkdir(savedir);  end

addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/Deeplearning_spike/Squiggles/func')

%%
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
% ask = [];
%
% clc
% disp('Select modality (MEG:1, EEG:2)')
% ask.modsel = input('');
%
% switch ask.modsel
%     case 1
%         modal = 'meg';
%     case 2
%         modal = 'eeg';
% end
%
% disp('Select frequncy range (Hz)');
% disp('1: Wideband 4-40');
% disp('2: Slow-rate 1-5')
% disp('3: High-rate 20-40');
% disp('4: Other frequncy');
% ask.freq_occur_sel = input(':');
%
% switch ask.freq_occur_sel
%     case 1
%         foi = [4,40];
%     case 2
%         foi = [1,5];
%     case 3
%         foi = [25,40];
%     case 4
%         disp('enter range [f1,f2]Hz:'); freq_range_sel = input(':');
%         foi = [freq_range_sel(1), freq_range_sel(2)];
% end


%%
pca_eeg = zeros(1,68);
pca_meg = zeros(1,68);
pca_all = zeros(1,68);

pca_eeg_sub = [];
pca_meg_sub = [];

ft_progress('init', 'etf', 'Please wait...');


for i= 1: length(d)
    disp([num2str(i),'/',num2str(length(d))])
    [pathstr, name] = fileparts(d(i).name);
    A = load(d(i).name);
    
    if isfield(A,'anot_data_all')
        
        %     ft_progress(i/length(d), 'Processing thresholded time windows %d from %d', i, length(d))  % show string, x=i/N
        %     pca_all = zeros(1,68);
        for j=1:length(A.anot_data_all)
            
            %         disp([num2str(j), '/', num2str(length(anot_data_all))])
            anot_data = A.anot_data_all{j};
            
            %         cfg = [];
            %         cfg.channel = 'EEG*';
            %         eeg = ft_selectdata(cfg, anot_data);
            %         %         figure, plot(smooth(mean(eeg.trial{:},1)))
            %
            %         cfg = [];
            %         cfg.channel = 'MEG*';
            %         meg = ft_selectdata(cfg, anot_data);
            %         figure, plot(smooth(mean(meg.trial{:},1)))
            
            %         avg_eeg(j,:) = zeros(1,68); D_eeg = smooth(mean(eeg.trial{:},1)); avg_eeg(j,1:length(D_eeg)) = D_eeg;
            %         avg_meg(j,:) = zeros(1,68); D_meg = smooth(mean(meg.trial{:},1)); avg_meg(j,1:length(D_meg)) = D_meg;
            
            %         pca_eeg(j,:) = zeros(1,68); D_eeg = smooth(smooth(do_pca(eeg.trial{:},1))); pca_eeg(j,1:length(D_eeg)) = D_eeg;
            %         pca_meg(j,:) = zeros(1,68); D_meg = smooth(smooth(do_pca(meg.trial{:},1))); pca_meg(j,1:length(D_meg)) = D_meg;
            %         pca_all(j,:) = zeros(1,68);
            D_all = smooth(smooth(do_pca(abs(anot_data.trial{:}),1)));
            pca_all(j,1:length(D_all)) = D_all;
            
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
        
        %     pca_eeg_sub {i} = pca_eeg;
        %     pca_meg_sub {i} = pca_meg;
        pca_all_sub {i} = pca_all;
    else
        disp('skipped')
    end
end

% figure,plot(mean((pca_all_sub {1}(:,1:67)),1))
% figure,plot(mean((pca_all_sub {2}(:,1:67)),1))
% figure,plot(mean((pca_all_sub {3}(:,1:67)),1))
% figure,plot(mean((pca_all_sub {4}(:,1:67)),1))
% figure,plot(mean((pca_all_sub {5}(:,1:67)),1))
% figure,plot(mean((pca_all_sub {6}(:,1:67)),1))

figure,plot(mean((pca_all(:,1:67)),1))

% figure,plot(mean((pca_meg(:,1:67)),1))
% % figure,plot(mean(abs(avg_meg(1:67)),1))

% pca_test = do_pca(pca_all(:,1:67),1);
% figure,plot(pca_test)
%%
size(pca_all)

data_pca = pca_all(:,1:67);

% pca_all_sub(1:5);
% numChannels = size(data{1},1);

figure
% tiledlayout(2,2)

for i = 1:20
    %     nexttile
    subplot(4,5,i)
    plot(data_pca(i,:))
    xlabel("Time Step")
end

numObservations = size(data_pca,1);
idxTrain = 1:floor(0.9*numObservations);
idxTest = floor(0.9*numObservations)+1:numObservations;
dataTrain = data_pca(idxTrain,:);
dataTest = data_pca(idxTest,:);

%%
XTrain = [];
for n = 1:size(dataTrain,1)
    X = dataTrain(n,:);
    XTrain{n} = X(:,1:end-1);
    TTrain{n} = X(:,2:end);
end

muX = mean(cat(2,XTrain{:}),2);
sigmaX = std(cat(2,XTrain{:}),0,2);

muT = mean(cat(2,TTrain{:}),2);
sigmaT = std(cat(2,TTrain{:}),0,2);

for n = 1:numel(XTrain)
    XTrain{n} = (XTrain{n} - muX) ./ sigmaX;
    TTrain{n} = (TTrain{n} - muT) ./ sigmaT;
end

layers = [
    sequenceInputLayer(numChannels)
    lstmLayer(128)
    fullyConnectedLayer(numChannels)
    regressionLayer];

options = trainingOptions('adam', 'MaxEpochs', 200, ...
    'Shuffle', 'every-epoch', ...
    'Plots','training-progress', ...
    'Verbose', 0);

% 'SequencePaddingDirection', 'left'

load('net')
% net = trainNetwork(XTrain,TTrain,layers,options);
% save('net','net')

%%
% link : https://www.mathworks.com/help/deeplearning/ug/time-series-forecasting-using-deep-learning.html
data = pca_all(:,1:67);
size(data)
numObservations = numel(data);
idxTrain = 1:floor(0.9*numObservations);
idxTest = floor(0.9*numObservations)+1:numObservations;
dataTrain = data(idxTrain);
dataTest = data(idxTest);

for n = 1:numel(dataTrain)
    X = dataTrain{n};
    XTrain{n} = X(:,1:end-1);
    TTrain{n} = X(:,2:end);
end


%%
layers = [ ...
    imageInputLayer([28 28 1])
    convolution2dLayer(5,20)
    reluLayer
    maxPooling2dLayer(2,'Stride',2)
    fullyConnectedLayer(10)
    softmaxLayer
    classificationLayer];
options = trainingOptions('sgdm', 'Plots', 'training-progress');
net = trainNetwork(XTrain, YTrain, layers, options);

%%
numTrainingFiles = 750;
[imdsTrain,imdsTest] = splitEachLabel(imds,numTrainingFiles,'randomize');

%%
% outsum_all = [];
% % for i= subid%1:length(d)
% for i= 1:20%length(d)
%     disp([num2str(i),'/',num2str(length(d))])
%     [pathstr, name] = fileparts(d(i).name);
%     load(d(i).name);
%
%     outsum = [];
%     for j=1:length(anot_data_all)
%         anot_data = anot_data_all{j};
%
%         cfg = [];
%         cfg.plot = 0;
%         outsum(j,:,:) = do_conn(cfg,anot_data.trial{1});
%
%         %         cfg = [];
%         %         cfg.blocksize = anot_data.time{1}(end) - anot_data.time{1}(1);
%         %         cfg.viewmode = 'vertical'; %butterfly';
%         %         cfg.continuous = 'yes';
%         %         cfg.axisfontsize = 7;
%         %         cfg.fontsize = 7;
%         %         cfg.channel = 'EEG*';
%         %         cfg.preproc.demean = 'yes';
%         %         cfg.position = [300   900   500   1500];
%         %         ft_databrowser(cfg, anot_data);
%         %         cfg.channel = 'MEG*';
%         %         cfg.position = [850   900   500   1500];
%         %         ft_databrowser(cfg, anot_data);
%         %
%         %         pause,
%         %         close all,
%     end
%     outsum_all{i} = squeeze(mean(outsum,1));
% end



%% Run training


% addpath('/data/MEG/Research/awang/Other tools/Kaggle-EEG')
% %
% params = [];
% % Set model and cv parameters
% % CV
% params.cvParams.cvMode = 'Custom';
% params.cvParams.k = 6;
% params.cvParams.evalProp = 0.2;
% params.cvParams.overSample = 0.05;
% params.cvParmas.seed = 2222;
% % Both models
% params.modParams.keepIdx = featuresTrain.keepIdx;
% params.modParams.prior = 'Empirical';
% params.modParams.hyper = 0;
% params.modParams.standardize = true;
% params.modParams.seed = 1111;
% % SVM
% params.modParams.polyOrder = 2;
% params.modParams.BC = 1000;
% % RBT
% params.modParams.nLearners = 100;
% params.modParams.LearnRate = 1;
% params.modParams.MaxNumSplits = 20;
%
% % Run train function
% [SVMg, RBTg] = trainModels(featuresTrain, params);

%%
%         anot_data
%         COEFS = cwt(s,1:32,'cgau4');
%         COEFS = cwt(anot_data.trial{1});
%
%
%         cfg = [];
%         cfg.output     = 'pow';
%         cfg.channel    = 'all';
%         cfg.method     = 'mtmconvol';
%         cfg.method     = 'wavelet';
%         cfg.taper      = 'hanning';
%         if foi(2) - foi(1) < 10
%             cfg.foi        = foi(1):1:foi(2);
%         else
%             cfg.foi        = foi(1):2:foi(2);
%         end
%         cfg.keeptrials = 'yes';
%         cfg.t_ftimwin  = 3 ./ cfg.foi;
%         cfg.tapsmofrq  = 0.8 * cfg.foi;
%         cfg.toi        = anot_data.time{1}(1):0.05:anot_data.time{1}(end);
%         tfr_data        = ft_freqanalysis(cfg, anot_data);
%
%         cfg = []; cfg.savepath = 1; cfg.savefile = [];
%         cfg.fmax = foi(2);
%         cfg.toi = [tfr_data.time(1), tfr_data.time(end)];
%         cfg.bslcorr = 2; cfg.plotflag = 1; cfg.title = modal;
%         [~,~, tfr_val]    = do_tfr_plot(cfg, tfr_data);