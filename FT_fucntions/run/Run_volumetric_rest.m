%% Volumetric-based analysis
% mridir = fullfile(indir,subj,'Anatomy/mri');
% d = rdir(fullfile(mridir,subj,'subjectimage*.mat'));

mridir = fullfile(indir,subj,'Anatomy/mri');
fidfile = fullfile(indir,subj,'Anatomy/bem/Anatomy-fiducials.fif');

% /rcc/stor1/projects/ECP/MEG/MEG_Work/EC1002/Anatomy/mri/T1.mgz
% /rcc/stor1/projects/ECP/MEG/MEG_Work/EC1002/Anatomy/bem/Anatomy-fiducials.fif

clear fid
% if ~isempty(d)
mripfile = fullfile(mridir,'T1.mgz');
if exist(outputmridir, 'file') == 0, mkdir(outputmridir); end

headshape = ft_read_headshape(datafile);

cfg = [];
cfg.megdata = ep_data.grad;
cfg.mripfile = mripfile;
cfg.hsfile = datafile; % headshape;
cfg.fid = headshape.fid;
cfg.outputmridir = outputmridir;
cfg.subj = subj;
cfg.plotflag = 2;
cfg.atlas = atlas;
%     [mri_realigned,individual_headmodel,headshape, individual_grid_8mm, individual_grid_10mm] = vy_mri_neuromag2(cfg);
[mri_realigned,individual_headmodel,headshape, individual_grid_8mm, individual_grid_10mm] = vy_mri_neuromag5(cfg);
%     vy_do_freesurfer(cfg);
% end
cd(outd.sub)


%% Choosing mesh
choose_grid = 2;
switch choose_grid
    % if flag.warping == 1
    case 1
        switch flag.meshgrid_sel
            case 1
                individual_grid = outanat.individual_grid_10mm_indiv;
            case 2
                individual_grid = outanat.individual_grid_8mm_indiv;
        end
    case 2
        switch flag.meshgrid_sel
            case 1
                meshtag = 'lowres';
                %         load('standard_sourcemodel3d10mm');
                load temp_grid % low-res
                template_grid = ft_convert_units(template_grid, 'mm');
                individual_grid = individual_grid_10mm;
            case 2
                meshtag = 'highres';
                %         load('standard_sourcemodel3d8mm');
                load temp_grid_8mm % high-res
                individual_grid = individual_grid_8mm;
        end
        % else
end

%% Anatomoy check!
saveflag = 2;
if flag.anatomy_check == 1
    cfg = [];
    cfg.saveflag = [];
    cfg.headmodel = individual_headmodel;
    cfg.leadfield = individual_grid;
    cfg.mri_realigned  = mri_realigned;
    cfg.headshape = headshape;
    cfg.outputmridir = outputmridir;
    cfg.mtd = 'vol';
    vy_mri_inspection(cfg, ep_data);
    %     vy_mri_inspection(t_data, individual_headmodel,individual_grid,headshape, mri_realigned,outputmridir,saveflag);
end

%%
template_mri = ft_read_mri(fullfile(allpath.ft_path,'template/anatomy','single_subj_T1.nii')); %

%%
cfg = [];
cfg.savefile = [];
cfg.saveflag = 2;
cfg.foilim = [2 50];
cfg.plotflag  = 2;
cfg.tapsmofrq       = 1;
cfg.taper    = 'hanning';
f_data = vy_fft(cfg, ep_data); f_data.elec = sens;

% % freq analysis - prepration for DICS source analysis
% f_data.bsl = vy_fft(ep_data.bsl, [2,40], 0,[],0); f_data.bsl.elec = sens;
% f_data.pst = vy_fft(ep_data.pst, [2,40], 0,[],0); f_data.pst.elec = sens;

% PSD - sensor space
psd = squeeze(mean(mean(abs(f_data.fourierspctrm),2),1));
ff = linspace(1, cfg.foilim(2), length(psd));

figure,plot(ff,psd)
xlabel('Hz'); ylabel('psd'),legend({'psd'})

%%
f = input('Freq of interest? ');
% f = 22;
% tapsmofrq = input('tapsmofrq? ');
tapsmofrq = 4;

cfg = [];
cfg.savefile = [];
cfg.saveflag = 2;
cfg.foilim = [f f];
cfg.plotflag  = 2;
cfg.taper    = 'dpss'; cfg.tapsmofrq  = tapsmofrq;

if f < 4, cfg.tapsmofrq  = 1; cfg.taper    = 'hanning'; end
if toi(1,2)-toi(1,1) < 0.4, cfg.taper    = 'hanning'; end

f_data = [];
[f_data.pst,~,~,~] = vy_fft(cfg, ep_data); f_data.pst.elec = datain.grad;
cfg.foilim = [15 15];
cfg.tapsmofrq  = 12;
[f_data.bsl,~,~,~] = vy_fft(cfg, ep_data); f_data.bsl.elec = datain.grad;

%                 cfg = [];
%                 f_data.app = ft_appenddata(cfg,f_data.bsl,f_data.pst);

%%
switch method
    case 1
        %% LCMV
        %         mtag = 'lcmv_stat';
        mtag = 'lcmv';
        outd.vol = fullfile(outd.sub, mtag);
        cfg = [];
        cfg.allpath = allpath;
        cfg.grid = individual_grid;
        cfg.headmodel = individual_headmodel;
        cfg.subj = subj;
        cfg.sens = sens;
        %         cfg.mtag = 'lcmv'; cfg.filterflag =  1;
        cfg.mtag = mtag; cfg.filterflag =  1;
        outd.vol = fullfile(outd.sub,cfg.mtag);
        %         cfg.fb   = [12,30]; % Hz
        %         cfg.fb   = [16,25]; % Hz
        cfg.fb   = [14,27]; % Hz
        cfg.atlas = allpath.atlas_path;
        cfg.outputdir = outd.vol;
        cfg.toi       = {toi};
        cfg.template_grid = template_grid;
        cfg.template_mri = template_mri;
        cfg.plotflag     = 1;
        switch mtag
            case 'lcmv'
                vy_source_lcmv(cfg, cln_data);
            case 'lcmv_stat'
                vy_source_lcmv_stats(cfg, cln_data);
        end
        
    case 2
        Run_rest_conn
        
    case 3
        %%
        %         disp('1: DICS, simple contrasting');
        %         disp('2: DICS-stats');
        %          mtag = input(':');
        mtag = 1;
        
        switch mtag
            case 1
                mtag = 'dics';
            case 2
                mtag = 'dics_stat';
        end
        
        %         mtag = 'dics_ratio';
        %         mtag = 'dics_fs';
        if rl == 1, mtag_lab = 'dics_res'; else, mtag_lab = 'dics'; end
        outd.vol = fullfile(outd.sub,mtag_lab);
        
        switch mtag
            case 'dics'
                
                
                %%
                %- FFT_based
                cfg = [];
                cfg.method = 'dics';
                cfg.dics.lambda = '100%';
                cfg.grid  = individual_grid;
                cfg.frequency    = f_data.pst.freq;
                cfg.headmodel = individual_headmodel;
                %                 cfg.dics.keepfilter = 'yes';
                cfg.dics.fixedori    = 'yes'; % project on axis of most variance using SVD
                %                 cfg.dics.projectnoise = 'yes';
                s_data.bsl      = ft_sourceanalysis(cfg, f_data.bsl);
                s_data.pst      = ft_sourceanalysis(cfg, f_data.pst);
                
                %                 cfg.dics.lambda       = 0;
                %                 sourceavg = ft_sourceanalysis(cfg, f_data.pst);
                
                %
                %                 sourceNAI = sourceavg;
                %                 sourceNAI.avg.pow = sourceavg.avg.pow ./ sourceavg.avg.noise;
                %                 source_diff_dics = sourceNAI;
                
                %% Nagative effects
                tmp = s_data.bsl.avg.pow; tmp(isnan(tmp))=0; tmp = tmp./max(tmp); s_data.bsl.avg.pow = tmp;
                tmp = s_data.pst.avg.pow; tmp(isnan(tmp))=0; tmp = tmp./max(tmp); s_data.pst.avg.pow = tmp;
                
                %                 cfg = [];
                %                 cfg.parameter = 'pow';
                %                 %     cfg.operation = '(x1-x2)/(x1+x2)';
                %                 cfg.operation = 'log10(x1/x2)'; % sourceA divided by sourceB
                %                 source_diff_dics = ft_math(cfg,s_data.pst,s_data.bsl);
                %                 source_diff_dics.pow(isnan(source_diff_dics.pow))=0;
                %                 source_diff_dics.pow = (s_data.pst.avg.pow - s_data.bsl.avg.pow)./(s_data.pst.avg.pow + s_data.bsl.avg.pow);
                % %                 source_diff_dics.pow(source_diff_dics.pow>0)=0;
                % %                 source_diff_dics.pow = abs(source_diff_dics.pow);
                %                 source_diff_dics.pos     = template_grid.pos;
                %                 source_diff_dics.dim     = template_grid.dim;
                %                 source_diff_dics.inside  = template_grid.inside;
                
                
                cfg = [];
                cfg.parameter = 'avg.pow';
                cfg.operation = '(x1-x2)';
                source_diff_dics = ft_math(cfg,s_data.pst,s_data.bsl);
                source_diff_dics.pow = (s_data.pst.avg.pow - s_data.bsl.avg.pow)./(s_data.pst.avg.pow + s_data.bsl.avg.pow);
                source_diff_dics.pow(source_diff_dics.pow<0)=0;
                
                cfg = [];
                cfg.mask = 'pow';
                %                 cfg.loc = 'max';
                cfg.template = mri_realigned;
                cfg.savefile = [];
                cfg.volnorm     = 2; % yes: 1
                cfg.method = 'ortho';
                source_int_all = vy_source_plot(cfg, source_diff_dics);
                
                %%
                disp('1: Yes');
                disp('2: No');
                sm = input('template + Surface map + laterality: ');
                
                if sm ==1
                    
                    cfg = [];
                    cfg.parameter = 'avg.pow';
                    cfg.operation = '(x1-x2)';
                    source_diff_dics_temp = ft_math(cfg,s_data.pst,s_data.bsl);
                    source_diff_dics_temp.pos     = template_grid.pos;
                    source_diff_dics_temp.dim     = template_grid.dim;
                    source_diff_dics_temp.inside  = template_grid.inside;
                    source_diff_dics_temp.pow = (s_data.pst.avg.pow - s_data.bsl.avg.pow)./(s_data.pst.avg.pow + s_data.bsl.avg.pow);
                    source_diff_dics_temp.pow(source_diff_dics_temp.pow<0)=0;
                    
                    %                 source_diff_dics.pow(source_diff_dics.pow>0)=0;
                    
                    cfg = [];
                    cfg.mask = 'pow';
                    cfg.loc = 'min';
                    cfg.template = template_mri;
                    cfg.savefile = 'dics1';
                    cfg.method = 'ortho';
                    cfg.volnorm     = 2; % yes: 1
                    source_int_dics_templ = vy_source_plot(cfg, source_diff_dics_temp);
                    
%                     savepath = [];
%                     savepath{1} = 'dics2.png';
%                     savepath{2} = 'dics3.png';
%                     
%                     cfg = [];
%                     cfg.subj = subj;
%                     cfg.mask = 'pow';
%                     cfg.thre = 0.5;
%                     cfg.savepath = savepath;
%                     cfg.colorbar = 2;
%                     vy_mapvisualisation(cfg, source_int_dics_templ);
%                     %                 view([90 90]);
%                     %                 camlight; material dull;
                    
                    mask = 'pow';
                    [~, D_par, coor] = vy_parcellate(source_diff_dics_temp, atlas, mask);
                    D_par.powdimord = 'chan';
                    [ROI, ROI_sel] = vy_ROI_report(D_par,.8, coor, mask);
                    
                    
                    idx1 = [];
                    idx1.central = [2,20];
                    idx1.frontal = 4:2:18;
                    idx1.subcor = [22:2:48,78];
                    idx1.Occ = 50:2:56;
                    idx1.pari= 58:2:70;
                    idx1.temp = 82:2:90;
                    
                    % idx = [idx_central,idx_frontal,idx_subcor,idx_Occ,idx_pari,idx_temp];
                    idx2         = [idx1.frontal,idx1.temp,idx1.pari];
                    %                 idx2         = [idx1.frontal, idx1.temp];
                    %                 idx2         = [idx1.frontal];
                    %                 idx2         = [idx1.temp];
                    
                    idx = [4:2:18,24:2:26,82:2:90];
                    
                    
                    data = [];
                    data.value = D_par.(mask);
                    data.label = D_par.label;
                    LI_meg = vy_laterality(data,idx);
                end
                
        end
end

