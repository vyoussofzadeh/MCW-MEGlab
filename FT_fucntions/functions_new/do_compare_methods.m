function do_compare_methods(cfg_main)

d_in = cfg_main.d_in;
network_sel = cfg_main.network_sel;

for LI_method=1:3
    
    switch LI_method
        case 1
            mlabel = 'threshold';
        case 2
            mlabel = 'counting';
        case 3
            mlabel = 'bootstrapping';
    end
    for contrast_type=1:4
        
        switch contrast_type
            case 1
                conlabel = 'LI_anim_hc';
            case 2
                conlabel = 'LI_anim_pt';
            case 3
                conlabel = 'LI_symb_hc';
            case 4
                conlabel = 'LI_symb_pt';
        end
        for j =1:length(network_sel)
            
            val = (d_in.(mlabel).(conlabel).LI_sub);
            
            (network_sel(j))
            mLI_sub1 = squeeze(mean(LI_anim_hc.LI_sub,2)); tag = [mlabel, '; anim hc'];
            
            clc, mLI_sub_hc = mLI_sub1;
            
            figure,
            clear LI_val
            for j=1:length(network_sel)
                hold on
                plot(mLI_sub_hc(network_sel(j),:),'LineWidth',3, 'color', colr(j,:)),
                val = round(mean(wi(:,1),2),2);
                set(gca,'Xtick', 1:2:length(wi),'XtickLabel',val(1:2:end));
                set(gca,'FontSize',8,'XTickLabelRotation',90);
                set(gcf, 'Position', [1000   400   1100   500]);
            end
            lgnd = legend([net_sel_mutiple_label(network_sel); 'mean']);
            title(tag)
            ylabel('LI')
            xlabel('time')
            set(gca,'color','none');
            set(lgnd,'color','none');
            
        end
    end
end
end

