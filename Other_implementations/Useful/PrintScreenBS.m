
% Orient = {'left'; 'right';'top';'bottom';'left_intern';'right_intern'};
% addpath('/data/MEG/Vahab/Github/MCW-MEGlab/Projects/ECP/Run')
% savedir = '/data/MEG/Vahab/test_data/test_resutls';
% FileName = 'Group_analysis/Beta/results_210420_1516.mat';
% cd('/group/jbinder/ECP/MEG/MEG_work_BS/data')
% a = load(FileName);
% savetag = a.Comment;
% Run_screenshot


%%

Orient = {'left'; 'right';'top';'bottom';'left_intern';'right_intern'};
% savedir = '/data/MEG/Research/spendl/spendl_SK_scripts/results';
savedir = '/data/MEG/Vahab/test_data/test_resutls';

% cd('/data/MEG/Research/spendl/Spendl_FTimport_SK_thetasync_250610/data')
cd('/group/jbinder/ECP/MEG/MEG_work_BS/data')


% FileName = 'Group_analysis/@intra/presults_fdr_210308_1726.mat';
FileName = 'Group_analysis/Beta/results_210420_1516.mat';

a = load(FileName);
% % Early_vs

% close all
% hFig = view_surface_data([], FileName, [], 'NewFigure');
% set(gca,'color','w');

close all
hFig = view_surface_data([], FileName, [], 'NewFigure');
set(gca,'color','w');
bst_colormaps('SetColorbarVisible', hFig, 0);
axis equal
pause,

for i=1:length(Orient)
    figure_3d('SetStandardView', hFig, Orient{i});
    img = out_figure_image(hFig, '', '');
    imgFile = fullfile(savedir, [Orient{i},'.tif']);
    out_image(imgFile, img);
    pause(1)
    %     saveas(gcf,imgFile)
    
end
disp(a.Comment)
disp('was saved at'),
disp(savedir)

%%
Orient = {'left'; 'bottom';'right';};

close all
hFig = view_surface_data([], FileName, [], 'NewFigure');
set(gca,'color','w');
bst_colormaps('SetColorbarVisible', hFig, 0);
axis equal

a = [];
for i=1:length(Orient)
    figure_3d('SetStandardView', hFig, Orient{i});
    img = out_figure_image(hFig, '', '');
%     imgFile = fullfile(savedir, [Orient{i},'.jpg']);
    imgFile = fullfile(savedir, 'all.tif');
    out_image(imgFile, img);
    b=imread(imgFile);
    if isa(b,'uint8'), b=double(b)/255; end
    if max(b(:))>1, b=double(b)/double(max(b(:))); end
    a{i}=double(b);
    %     waitbar(n1/numel(mosaic_commands),hw);
    pause(1)
    %     saveas(gcf,imgFile)
    
end

ncut = 16;
domosaiccrop = 1;
if domosaiccrop
    cropt_idx={};
    for n=1:numel(a)
        cropt=any(any(diff(a{n},1,2),2),3);
        cropt_idx{n,1}=max(1,sum(~cumsum(cropt))-ncut):size(a{n},1)-max(0,sum(~cumsum(flipud(cropt)))-ncut);
        cropt=any(any(diff(a{n},1,1),1),3);
        cropt_idx{n,2}=max(1,sum(~cumsum(cropt))-ncut):size(a{n},2)-max(0,sum(~cumsum(flipud(cropt)))-ncut);
    end
end


if domosaiccrop
    cropt_idx1234=cropt_idx{1,1}; for n=2:numel(a), cropt_idx1234=union(cropt_idx1234,cropt_idx{n,1}); end
    ta=[]; for n=1:numel(a), ta=cat(2,ta,a{n}(cropt_idx1234,cropt_idx{n,2},:)); end; a=ta;
    %cropt_idx1234=union(union(union(cropt_idx{1,1},cropt_idx{2,1}),cropt_idx{3,1}),cropt_idx{4,1});
    %a=[a{1}(cropt_idx1234,cropt_idx{1,2},:),a{2}(cropt_idx1234,cropt_idx{2,2},:),a{3}(cropt_idx1234,cropt_idx{3,2},:),a{4}(cropt_idx1234,cropt_idx{4,2},:)];
else
    a=cat(2,a{:});
    %a=[a{1},a{2},a{3},a{4}];
end

imwrite(a,imgFile);
% % delete(hw);
