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
        %%
        cfg = [];
        mtag = 'conn';
        %         mtag = 'conn_bl'; cfg.fl = [1, 10];% band limited
        %         mtag = 'conn_bs'; cfg.fb = 10; % band-stop
        outd.vol = fullfile(outd.sub,mtag);
        cfg.allpath = allpath;
        cfg.grid = individual_grid;
        cfg.headmodel = individual_headmodel;
        cfg.subj = subj;
        cfg.sens = sens;
        cfg.mtag = mtag;
        cfg.toi  = toi;
        cfg.atlas = allpath.atlas_path;
        cfg.outputdir = outd.vol;
        cfg.template_grid = template_grid;
        cfg.template_mri = template_mri;
        switch mtag
            case 'conn'
                vy_network_light1(cfg,ep_data) % conn-network analysis
            case {'conn_bs','conn_bl'}
                vy_network_light1(cfg, cln_data) % conn-network analysis
        end
        
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
                
                %% Revert to new ft!
%                 restoredefaultpath
%                 addpath((allpath.ft_path));
%                 ft_defaults
%                 addpath(genpath(allpath.hcp_path));
%                 addpath(genpath(allpath.cd_org));
%                 addpath(genpath(allpath.exfig_path));
                
                %%
                f = input('Freq of interest? ');
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
                %- FFT_based
%                 cfg = [];
%                 cfg.method = 'dics';
%                 cfg.dics.lambda = '100%';
%                 cfg.grid  = individual_grid;
%                 cfg.frequency    = f_data.pst.freq;
%                 cfg.headmodel = individual_headmodel;
%                 %                 cfg.dics.keepfilter = 'yes';
%                 cfg.dics.fixedori    = 'yes'; % project on axis of most variance using SVD
% %                 cfg.dics.projectnoise = 'yes';
%                 s_data.bsl      = ft_sourceanalysis(cfg, f_data.bsl);
%                 s_data.pst      = ft_sourceanalysis(cfg, f_data.pst);
                
                %                 cfg.dics.lambda       = 0;
                %                 sourceavg = ft_sourceanalysis(cfg, f_data.pst);
                
                %
                %                 sourceNAI = sourceavg;
                %                 sourceNAI.avg.pow = sourceavg.avg.pow ./ sourceavg.avg.noise;
                %                 source_diff_dics = sourceNAI;
                
                %%
%                 restoredefaultpath
%                 addpath(genpath(allpath.ft_old));
%                 addpath(genpath(allpath.hcp_path));
%                 addpath(genpath(allpath.cd_org));
%                 ft_defaults
                
                %%
                %- FFT_based
                cfg = [];
                cfg.method = 'pcc';
                cfg.pcc.lambda = '5%';
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
                %                 gtm = 'degrees';
                
                figure,imagesc(conn_pst.(mtd) - conn_bsl.(mtd)); colorbar
%                 figure,imagesc(conn_bsl.(mtd)); colorbar
%                 figure,imagesc(conn_pst.(mtd)); colorbar
                
                %                 conn_bsl.(mtd)(isnan(conn_bsl.(mtd))) = 0;
                %                 conn_pst.(mtd)(isnan(conn_pst.(mtd))) = 0;
                %
                %                 tmp = conn_bsl.(mtd); tmp(isnan(tmp(:)))=0; tmp = tmp./max(tmp(:)); conn_bsl.(mtd) = tmp;
                %                 tmp = s_data.pst.avg.pow; tmp(isnan(tmp))=0; tmp = tmp./max(tmp(:)); s_data.pst.avg.pow = tmp;
                
                
                tmp = conn_bsl.(mtd); tmp = tmp./max(tmp(:)); conn_bsl.(mtd) = tmp;
                tmp = conn_pst.(mtd); tmp = tmp./max(tmp(:)); conn_pst.(mtd) = tmp;
                
                cfg = [];
                cfg.method    = gtm;
                cfg.parameter = mtd;
                cfg.threshold = 0.7.*max(conn_bsl.(mtd)(:));
                net_bsl = ft_networkanalysis(cfg,conn_bsl);
                net_bsl.pos = s_data.bsl.pos;
                net_bsl.dim  = s_data.bsl.dim;
                net_bsl.inside  = s_data.bsl.inside;
                
                cfg = [];
                cfg.method    = gtm;
                cfg.parameter = mtd;
                cfg.threshold = 0.7.*max(conn_pst.(mtd)(:));
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
                
                cfg = [];
                cfg.mask = gtm;
                cfg.loc = 'min';
                cfg.template = template_mri;
                cfg.savefile = 'conn1';
                cfg.volnorm     = 2; % yes: 1
                net_diff1 = vy_source_plot(cfg, net_diff);
                
                %%
                savepath = [];
                savepath{1} = 'conn2.png';
                savepath{2} = 'conn3.png';
                
                cfg = [];
                cfg.subj = subj;
                cfg.mask = gtm;
                cfg.thre = 0.6;
                cfg.savepath = savepath;
                vy_mapvisualisation(cfg, net_diff1);
                
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
                source_diff_dics.pos     = template_grid.pos;
                source_diff_dics.dim     = template_grid.dim;
                source_diff_dics.inside  = template_grid.inside;
                source_diff_dics.pow = (s_data.pst.avg.pow - s_data.bsl.avg.pow)./(s_data.pst.avg.pow + s_data.bsl.avg.pow);
                source_diff_dics.pow(source_diff_dics.pow<0)=0;

                %                 source_diff_dics.pow(source_diff_dics.pow>0)=0;                
                
                cfg = [];
                cfg.mask = 'pow';
                cfg.loc = 'min';
                cfg.template = template_mri;
                cfg.savefile = 'dics1';
                cfg.volnorm     = 2; % yes: 1
                source_int_dics = vy_source_plot(cfg, source_diff_dics);
                
                savepath = [];
                savepath{1} = 'dics2.png';
                savepath{2} = 'dics3.png';
                
                cfg = [];
                cfg.subj = subj;
                cfg.mask = 'pow';
                cfg.thre = 0.6;
                cfg.savepath = savepath;
                vy_mapvisualisation(cfg, source_int_dics);
%                 view([90 90]); 
%                 camlight; material dull;

        end
end

