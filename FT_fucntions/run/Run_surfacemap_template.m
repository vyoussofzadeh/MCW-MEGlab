% disp('1: Yes')
% disp('2: No')
% surask = input('Surface map:');
% switch surask
%     case 1
        s_in.pos     = template_grid.pos;
        s_in.dim     = template_grid.dim;
        s_in.inside  = template_grid.inside;
        
        cfg = [];
        cfg.mask = mask;
        cfg.loc = 'max';
        cfg.template = template_mri;
        cfg.savefile = saveID;
        cfg.volnorm     = 2; % yes: 1
        cfg.method        = 'ortho';
        cfg.nslices = [10,60];
        source_int = vy_source_plot(cfg, s_in);
        
        cfg = [];
        cfg.subj = saveID;
        cfg.mask = mask;
        cfg.thre = thre;
        cfg.savepath = 1;
        cfg.colorbar = 2;
        vy_mapvisualisation(cfg, source_int);
        % vy_mapvisualisation(source_int_dics,cfg.mask,0.6, []);
        
% end