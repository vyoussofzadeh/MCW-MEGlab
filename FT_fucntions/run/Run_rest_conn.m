%- FFT_based
cfg = [];
cfg.method = 'pcc';
cfg.pcc.lambda = '100%';
cfg.senstype  = 'MEG';
cfg.grid  = individual_grid;
cfg.frequency = f_data.pst.freq;
cfg.headmodel = individual_headmodel;
%                 cfg.dics.keepfilter = 'yes';
cfg.pcc.fixedori = 'yes'; % project on axis of most variance using SVD
%                 cfg.dics.projectnoise = 'yes';
s_data.bsl      = ft_sourceanalysis(cfg, f_data.bsl);
s_data.pst      = ft_sourceanalysis(cfg, f_data.pst);

%%
disp('1:coh');
disp('2:imgcoh');
disp('3:plv');
disp('4:powcorr');
disp('5:powcorr_ortho');
disp('6:wpli_debiased');
conn_sel = input(':');
cfg           = [];
switch conn_sel
    case 1
        cfg.method    = 'coh'; cfg.complex = 'abs'; label = 'coh'; savetag = 'coh';
    case 2
        cfg.method    = 'coh'; cfg.complex = 'absimag'; label = 'coh'; savetag = 'imcoh';
    case 3
        cfg.method    = 'plv'; cfg.complex = 'abs'; label = 'plv'; savetag = 'plv';
    case 4
        cfg.method    = 'powcorr'; cfg.complex = 'abs'; label = 'powcorr'; savetag = 'powcorr';
    case 5
        cfg.method    = 'powcorr_ortho'; cfg.complex = 'abs'; label = 'powcorr'; savetag = 'powcorr-ortho';
    case 6
        cfg.method    = 'wpli_debiased'; cfg.complex = 'abs'; label = 'wpli_debiased'; savetag = 'wpli-debiased';
end
% cfg.method    = 'wpli'; cfg.complex = 'abs'; label = 'wpli';
% cfg.method    = 'psi'; cfg.complex = 'abs'; label = 'psi'; cfg.bandwidth = 4;
% cfg.method    = 'wppc'; cfg.complex = 'abs'; label = 'wppc';


tmp = cat(1, s_data.bsl.avg.mom{s_data.bsl.inside});
f_data1.bsl = f_data.bsl;
f_data1.bsl.fourierspctrm = tmp';
lbl = [];
for i=1:size(f_data1.bsl.fourierspctrm,2)
    lbl{i}=num2str(i);
end
f_data1.bsl.label = lbl';


tmp = cat(1, s_data.pst.avg.mom{s_data.pst.inside});
f_data1.pst = f_data.pst;
f_data1.pst.fourierspctrm = tmp';
lbl = [];
for i=1:size(f_data1.pst.fourierspctrm,2)
    lbl{i}=num2str(i);
end
f_data1.pst.label = lbl';

conn_bsl           = ft_connectivityanalysis(cfg, f_data1.bsl); conn_bsl.dimord    = 'pos_pos';
conn_pst           = ft_connectivityanalysis(cfg, f_data1.pst); conn_pst.dimord    = 'pos_pos';

%                 conn_bsl           = ft_connectivityanalysis(cfg, s_data.bsl); conn_bsl.dimord    = 'pos_pos';
%                 conn_pst           = ft_connectivityanalysis(cfg, s_data.pst); conn_pst.dimord    = 'pos_pos';

%%
mtd = [label,'spctrm']; % method par.
gtm = 'eigenvector_cent'; % grapth theory measure
% gtm = 'degrees';

% figure,imagesc(conn_pst.(mtd) - conn_bsl.(mtd)); colorbar
%                 figure,imagesc(conn_bsl.(mtd)); colorbar
%                 figure,imagesc(conn_pst.(mtd)); colorbar

%                 conn_bsl.(mtd)(isnan(conn_bsl.(mtd))) = 0;
%                 conn_pst.(mtd)(isnan(conn_pst.(mtd))) = 0;
%
%                 tmp = conn_bsl.(mtd); tmp(isnan(tmp(:)))=0; tmp = tmp./max(tmp(:)); conn_bsl.(mtd) = tmp;
%                 tmp = s_data.pst.avg.pow; tmp(isnan(tmp))=0; tmp = tmp./max(tmp(:)); s_data.pst.avg.pow = tmp;


tmp = conn_bsl.(mtd); tmp = tmp./max(tmp(:)); conn_bsl.(mtd) = tmp;
tmp = conn_pst.(mtd); tmp = tmp./max(tmp(:)); conn_pst.(mtd) = tmp;

thre = 0.5;

cfg = [];
cfg.method    = gtm;
cfg.parameter = mtd;
cfg.threshold = thre.*max(conn_bsl.(mtd)(:));
net_bsl = ft_networkanalysis(cfg,conn_bsl);
net_bsl.pos = s_data.bsl.pos;
net_bsl.dim  = s_data.bsl.dim;
net_bsl.inside  = s_data.bsl.inside;

cfg = [];
cfg.method    = gtm;
cfg.parameter = mtd;
cfg.threshold = thre.*max(conn_pst.(mtd)(:));
net_pst = ft_networkanalysis(cfg,conn_pst);
net_pst.pos = s_data.pst.pos;
net_pst.dim  = s_data.pst.dim;
net_pst.inside  = s_data.pst.inside;

%%
cfg = [];
cfg.parameter = gtm;
cfg.operation = 'x1-x2';
net_diff = ft_math(cfg,net_pst,net_bsl);

net_diff.(gtm) = zeros(size(net_diff.inside,1),1);
net_diff.(gtm)(net_diff.inside) = net_pst.(gtm) - net_bsl.(gtm);

%                 net_diff.(gtm) = net_pst.(gtm) - net_bsl.(gtm);
% net_diff.(gtm)(net_diff.(gtm)<0)=0;
% net_diff.(gtm)(net_diff.(gtm)>0)=0;
% net_diff = net_pst;
%                 net_diff.elec = freq_sel1.elec;

net_diff.pos     = template_grid.pos;
net_diff.dim     = template_grid.dim;
net_diff.inside  = template_grid.inside;
%                 net_diff.pow = (s_data.pst.avg.pow - s_data.bsl.avg.pow)./(s_data.pst.avg.pow + s_data.bsl.avg.pow);
net_diff.(gtm)(net_diff.(gtm)<0)=0;

%                 source_diff_dics.pow(source_diff_dics.pow>0)=0;

% savefig = savetag);

cfg = [];
cfg.mask = gtm;
cfg.loc = 'min';
cfg.template = template_mri;
cfg.savefile = savetag;
cfg.volnorm     = 2; % yes: 1
net_diff1 = vy_source_plot(cfg, net_diff);

%%
% savepath = [];
% savepath{1} = 'conn2.png';
% savepath{2} = 'conn3.png';

%%
clear savepath
savepath{1} = fullfile(['R_',savetag]);
savepath{2} = fullfile(['L_',savetag]);

% cfg = [];
% cfg.subj = subj;
% cfg.mask = gtm;
% cfg.thre = 0.6;
% cfg.savepath = savepath;
% cfg.colorbar = 2;
% vy_mapvisualisation(cfg, net_diff1);

