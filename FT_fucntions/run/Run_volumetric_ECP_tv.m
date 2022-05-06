

         mtag = 'dics';
         %         mtag = 'dics_ratio';
         %                 mtag = 'dics_stat';
         %         mtag = 'dics_fs';
         if rl == 1, mtag_lab = 'dics_res'; else, mtag_lab = 'dics'; end
         outd.vol = fullfile(outd.sub,mtag_lab);
         
         switch mtag
             case 'dics'
                 cfg = [];
                 cfg.grid = individual_grid;
                 cfg.allpath = allpath;
                 cfg.freq_of_interest  = freq_of_interest; % Hz
                 cfg.headmodel = individual_headmodel;
                 cfg.sens = sens;
                 cfg.mtag = mtag;
                 cfg.subj = subj;
                 cfg.toi = toi;
                 cfg.outputdir = outd.vol;
                 switch choose_grid
                     case 1
                         cfg.template_grid = [];
                         cfg.mri_aligned = outanat.mri_realigned;
                     case 2
                         cfg.template_grid = template_grid;
                 end
                 cfg.template_mri = template_mri;
                 cfg.fmax = fmax;
                 cfg.savedata = fullfile(outd.vol,[mtag,'_',subj]);
                 cfg.flag = flag;
                 cfg.plotting = 1;
                 vy_source_dics(cfg, ep_data);
                 
             case 'dics_ratio'
                 
                mtag = 'dics_ratio';
                cfg = [];
                cfg.grid = individual_grid;
                cfg.allpath = allpath;
                cfg.freq_of_interest   = freq_of_interest; % Hz
                cfg.headmodel = individual_headmodel;
                cfg.sens = sens;
                cfg.mtag = mtag;
                cfg.subj = subj;
                cfg.toi = toi;
                cfg.outputdir = outd.vol;
                cfg.template_grid = template_grid;
                cfg.template_mri = template_mri;
                cfg.savedata = fullfile(outd.vol,[mtag,'_',subj,'.mat']);
                vy_source_dics_ratio(cfg, ep_data);
                
            case 'dics_fs'
                %%
                anatomy_dir     = '/data/MEG/Clinical/ft_process/19/bednar_peggy/anat/BAK';
                load(fullfile(anatomy_dir,[subj,'_headmodel.mat']));
                load(fullfile(anatomy_dir,[subj,'_leadfield.mat']));
                load(fullfile(anatomy_dir,[subj,'_sourcemodel.mat']));
                
                if anatomy_check_flag == 1
                    cfg1 = [];
                    cfg1.saveflag = 2;
                    cfg1.headmodel = headmodel;
                    cfg1.sourcemodel = sourcemodel;
                    cfg1.leadfield = leadfield;
                    cfg1.mri = anatomy_dir;
                    cfg1.mtd = 'surf';
                    cfg1.headshape = headshape;
                    vy_mri_inspection(cfg1, t_data);
                    %                 individual_headmodel,individual_grid,headshape, mri_realigned,outputmridir,saveflag
                end
                cfg = [];
                cfg.headmodel = headmodel;
                cfg.sourcemodel = sourcemodel;
                cfg.leadfield   = leadfield;
                cfg.mtag = mtag;
                cfg.sens = sens;
                cfg.subj = subj;
                cfg.outputdir = outd.vol;
                vy_source_dics_fs(cfg, ep_data);
                
            case 'dics_stat'
                mtag = 'dics_stat';
                cfg = [];
                cfg.grid = individual_grid;
                cfg.allpath = allpath;
%                 cfg.f   = 20; % Hz
                cfg.headmodel = individual_headmodel;
                cfg.sens = sens;
                cfg.mtag = mtag;
                cfg.subj = subj;
                cfg.toi = toi;
                cfg.outputdir = outd.vol;
                cfg.template_grid = template_grid;
                cfg.template_mri = template_mri;
                vy_source_dics_stats(cfg, ep_data);
        end
        
end
