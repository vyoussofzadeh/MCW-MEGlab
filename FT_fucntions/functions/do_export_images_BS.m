function do_export_images_BS(cfg_in)


% clc, close all,

%%
disp('Where the map is goingn to be exported to ...')
% svdir = input('1: enter saving dir:');svdir = strrep(svdir, ' ', '');
svdir = cfg_in.svdir;
disp('===========')


disp('righ-click on the source map file/File/view file history')
% BSpath = input('2: enter BS path:');
BSpath = cfg_in.BSpath;
BSpath = strrep(BSpath, ' ', '');

% fname = input('3: enter BS source FileName:');
fname = cfg_in.fname;
fname = strrep(fname, ' ', '');
disp('===========')

disp('4: Adjust/check the threshold, press enter to proceed ...')
cd(BSpath)


%%
svname = cfg_in.sname;
if isempty(svname)
    sfile = load(fname); svname = sfile.Comment;
    disp(['suggesting name:',svname]);
    name_sel = input('1-suggested name, 2-other names:');
    if name_sel == 2
        svname = input('enter saving name:','s');
    else
        svname = sfile.Comment;
    end
end

%%
clc
disp('1: left;right;top;bottom;left_intern;right_intern')
disp('2: bottom;top;left;right;left_intern;right_intern')
disp('3: left;right;top;bottom;right_intern;left_intern')
disp('4: left;bottom;right')
disp('5: left;right;top;bottom')
disp('6: optional, e.g, {left;right;top}')
disp('7: left;right')
% side_sel = input(':');
side_sel =  cfg_in.side_sel;

switch side_sel
    case 1
        Orient = {'left'; 'right';'top';'bottom';'left_intern';'right_intern'};
    case 2
        Orient = {'bottom';'top';'left';'right';'left_intern';'right_intern'};
    case 3
        Orient = {'left'; 'right';'top';'bottom';'right_intern';'left_intern'};
    case 4
        Orient = {'left'; 'bottom';'right'};
    case 5
        Orient = {'left'; 'right';'top';'bottom'};
    case 6
        side_sel_man = input('enter selected views, in quotation marks');
        Orient = side_sel_man;
    case 7
        Orient = {'left';'right'};
end

%%
close all
hFig = view_surface_data([], fname, [], 'NewFigure');
set(hFig,'color','w');
bst_colormaps('SetColorbarVisible', hFig, 0);
axis equal

disp(cfg_in.seltime)
panel_time('SetCurrentTime', cfg_in.seltime);

% pause,

printoptions={'-djpeg90','-r600','-opengl'}; imgFile = fullfile(svdir, [svname,'.jpg']);
% printoptions={'-dpng', '-r300'}; imgFile = fullfile(svdir, [svname,'.png']);

a = [];
for i=1:length(Orient)
    figure_3d('SetStandardView', hFig, Orient{i});
    bst_colormaps('SetColorbarVisible', hFig, 0);
    
    switch cfg_in.imgres
        case 1
            %% high-res (Added by VYZ, 08/28/22)
            drawnow; print(hFig,printoptions{:},imgFile);
            b = imread(imgFile);
            
        case 2
            %% low-res print
            %     imgFile = fullfile(svdir, [svname,'.tif']);
            img = out_figure_image(hFig, '', '');
            imgFile = fullfile(svdir, [svname,'.tif']);
            out_image(imgFile, img);
            b = imread(imgFile);
    end
%%
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
disp('5: completed!, images were saved at,')
disp(svdir)

%%
end
