function do_peakcoor_detection(cfg)


%%
disp('Where the map is goingn to be exported to ...')
svdir = input('1: enter saving dir:','s');svdir = strrep(svdir, ' ', '');
disp('===========')


disp('righ-click on the source map file/File/view file history')
BSpath = input('2: enter BS path:','s'); BSpath = strrep(BSpath, ' ', '');
fname = input('3: enter BS source FileName:','s'); fname = strrep(fname, ' ', '');
disp('===========')

disp('4: Adjust/check the threshold, press enter to proceed ...')
cd(BSpath)

%%
sfile = load(fname);
svname = sfile.Comment;
disp(['suggesting name:',svname]);
svname = input('enter saving name:','s');
disp('===========')

%%
close all
npeaks = cfg.npeaks;

%%
cd(BSpath)
sSurface = load(fname);
cd ..
cd ./anat
d = load(sSurface.SurfaceFile);

%%
ImageGridAmp = sSurface.ImageGridAmp;
[maxval, maxidx] = max(ImageGridAmp);

b = [];
b.Vertices = d.Vertices;
b.Faces = d.Faces;
% figure; ft_plot_mesh(b, 'edgecolor', 'none', 'facealpha', 0.4);
% view([-180,0])
% hold on
% plot3(b.Vertices(maxidx,1),b.Vertices(maxidx,2),b.Vertices(maxidx,3),'*')

%%
pcoor = b.Vertices(maxidx,:);
coor = 1000.*b.Vertices;
[sor_val, sor_idx] = sort(ImageGridAmp,'descend');

%%
seed_coor = coor(sor_idx(1),:);

hFig = view_surface_data([], fname, [], 'NewFigure');
set(hFig,'color','w');
bst_colormaps('SetColorbarVisible', hFig, 0);
axis equal
alpha(0.5)

seed_val = [];
marker_size = 80;
hold on
x = b.Vertices(sor_idx(1),1); y = b.Vertices(sor_idx(1),2); z = b.Vertices(sor_idx(1),3);
kk = 1;
clear fh
fh(kk) = plot3(x,y,z,'m.','MarkerSize',marker_size);

r1 = x;
str = sprintf('%s', num2str(size(seed_coor,1)));

h(kk)   = text(x, y, z, str, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Interpreter', 'none','FontSize',10);
seed_val(kk) = ImageGridAmp(sor_idx(1));

for i=1:length(ImageGridAmp)-1
    if size(seed_coor,1) < npeaks
        d2 = [];
        for j=1:size(seed_coor,1)
            d2(j) = pdist2(seed_coor(j,:),coor(sor_idx(i+1),:));
        end
        md2 = min(d2);
        if md2 > 30
            seed_coor = [seed_coor;coor(sor_idx(i+1),:)];
            x = b.Vertices(sor_idx(i+1),1); y = b.Vertices(sor_idx(i+1),2); z = b.Vertices(sor_idx(i+1),3);
            kk = kk+1;
            seed_val(kk) = ImageGridAmp(sor_idx(i+1));
            fh(kk) = plot3(x,y,z,'m.','MarkerSize',marker_size);
            str = sprintf('%s', num2str(size(seed_coor,1)));
            h(kk)   = text(x, y, z, str, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Interpreter', 'none', 'FontSize',10);
        end
    end
end

%% Coordniates
sMRI = load('./@default_subject/subjectimage_T1.mat');

P_mni = [];
for k=1:length(seed_coor)
    P_mni(k,:) = 1e3*cs_convert(sMRI, 'scs', 'mni', seed_coor(k,:)/1e3);
end

disp('mni coordinates')
disp(round(P_mni))
rP_mni = round(P_mni);
disp('value')
disp(seed_val');

%%
T1 = table(round(P_mni));
T1.Properties.VariableNames{'Var1'} = 'mni';
T2 = table(round(seed_val'*100)/100);
T2.Properties.VariableNames{'Var1'} = 'value';
T = [T2,T1]; clear T1 T2

disp(fullfile(svdir))
disp('as,')
disp([svname,'.tif'])

writetable(T,[fullfile(svdir,svname),'.txt']);

%% left and right coordinates,
idx_left = find(rP_mni(:,1)<0);
idx_right = find(rP_mni(:,1)>0);

%% update color,
for j=1:length(idx_left), set(fh(idx_left(j)),'Color',[.98 .45 .02]),end
for j=1:length(idx_right), set(fh(idx_right(j)),'Color',[.88 .15 .42]),end

%% Apply trasnparency
for j=1:length(idx_left), alpha(fh(kk),0.5) ,end
for j=1:length(idx_right), alpha(fh(kk),0.5) ,end

%% 
disp('==================================================')
disp('1: left;right;top;bottom;left_intern;right_intern')
disp('2: bottom;top;left;right;left_intern;right_intern')
disp('3: left;right;top;;left_intern;right_intern')
disp('4: left;bottom;right')
disp('5: left;right;top;bottom')
disp('6: optional, e.g, {left;right;top}')
side_sel = input(':');

switch side_sel
    case 1
        Orient = {'left'; 'right';'top';'bottom';'left_intern';'right_intern'};
    case 2
        Orient = {'bottom';'top';'left';'right';'left_intern';'right_intern'};
    case 3
        Orient = {'left'; 'right';'top';'left_intern';'right_intern'};
    case 4
        Orient = {'left'; 'bottom';'right'};
    case 5
        Orient = {'left'; 'right';'top';'bottom'};
    case 6
        side_sel_man = input('enter selected views,');
        Orient = side_sel_man;
end

%% Export maps
alphalevel = input('sel transparecny 0-1(low-to-high): ');
alpha(alphalevel)

a = [];
for i=1:length(Orient)
    figure_3d('SetStandardView', hFig, Orient{i});
    panel_surface('SetSurfaceTransparency', hFig, 1, alphalevel);
    img = out_figure_image(hFig, '', '');
    imgFile = fullfile(svdir, [svname,'.tif']);
    out_image(imgFile, img);
    b=imread(imgFile);
    if isa(b,'uint8'), b=double(b)/255; end
    if max(b(:))>1, b=double(b)/double(max(b(:))); end
    a{i}=double(b);
    pause(1)
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
else
    a=cat(2,a{:});
end

imwrite(a,imgFile);
disp('5: completed!, images were saved at,')
disp(fullfile(svdir))
disp('as,')
disp([svname,'.tif'])
