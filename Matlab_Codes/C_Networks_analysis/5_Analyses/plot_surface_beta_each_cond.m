
addpath('C:\Users\Ness\Documents\Toolboxes\fieldtrip-20170116')
ft_defaults

load('E:\data_SELFI\last_analyses\icasso_nica15_corr25\pval_plot.mat');
load('E:\data_SELFI\last_analyses\icasso_nica15_corr25\betas_cond1_cond2.mat')

template_grid = load(['C:\Users\Ness\Documents\Toolboxes\fieldtrip-20170116\template\sourcemodel\standard_sourcemodel3d10mm.mat']);
templatedir  = ['C:\Users\Ness\Documents\Toolboxes\fieldtrip-20170116\external\spm8\templates'];
template = ft_read_mri(fullfile(templatedir, 'T1.nii'));
template_resliced = ft_volumereslice([], template);

%Max and min values
max_cond1 = max(mean(betas.cond1(16, :, :), 3))
min_cond1 = min(mean(betas.cond1(16, :, :), 3))

max_cond2 = max(mean(betas.cond2(16, :, :), 3))
min_cond2 = min(mean(betas.cond2(16, :, :), 3))

if max_cond1 > max_cond2
    max_all = max_cond1
else max_all = max_cond2
end

if min_cond1 < min_cond2
    min_all = min_cond1
else min_all = min_cond2
end

%Condition to plot
cond = 'Low contact'
comp15_betas.cond = betas.cond2(16, :, :)

plot_beta = mean(comp15_betas.cond, 3);

% Create a fieldtrip-like data structure with the betas values
Beta_cond_plot = [];
Beta_cond_plot.dim = template_grid.sourcemodel.dim;
Beta_cond_plot.inside = template_grid.sourcemodel.inside;
Beta_cond_plot.pos = template_grid.sourcemodel.pos;
Beta_cond_plot.method = 'average';
Beta_cond_plot.time{1} = pval_plot.time %data_hilbert_all.time;
Beta_cond_plot.avg.pow = plot_beta; % here are the p-values for each grid point
            
% Interpolation: Project on standard MRI brain
cfg             = [];
cfg.parameter   = 'pow';
Beta_cond_plot_intp  = ft_sourceinterpolate(cfg, Beta_cond_plot, template_resliced);

           % Left hemisphere internal side
            cfg              = [];
            cfg.method       = 'surface';
            cfg.funcolorlim  = [min_all abs(min_all)] %'maxabs' %'minzero'%'maxabs'
            cfg.funparameter = 'pow';
            cfg.camlight     = 'yes';
            cfg.colorbar = 'no'
            cfg.funcolormap   =  'jet';            
            cfg.surffile     = ['C:\Users\Ness\Documents\Toolboxes\fieldtrip-20170116\template\anatomy\surface_white_left.mat']
            ft_sourceplot(cfg, Beta_cond_plot_intp);
            view ([90 0]);
            set(gcf,'color','black')
            title(['Left hemisphere ' cond],'Color', 'w')
%             set(gca, 'CLim', [2.1 4.8])
%             colormap(mymap_GBOY)
            hc = colorbar;
            set(hc, 'Color', [1 1 1], 'FontName', 'cambria', 'box', 'off','Location', 'southoutside', 'FontWeight', 'bold')   %color bar in white
                         
            % Left hemisphere external side
            cfg              = [];
            cfg.method       = 'surface';
            cfg.funcolorlim  = [min_all abs(min_all)]
            cfg.funparameter = 'pow';
            cfg.camlight     = 'yes';
            cfg.colorbar = 'no'
            cfg.funcolormap   =  'jet';          
            cfg.surffile     = ['C:\Users\Ness\Documents\Toolboxes\fieldtrip-20170116\template\anatomy\surface_white_left.mat']
            ft_sourceplot(cfg, Beta_cond_plot_intp);
            view ([-90 0]);
            set(gcf,'color','black')
            title(['Left hemisphere ' cond],'Color', 'w')
%             set(gca, 'CLim', [2.1 4.8])
%             colormap(mymap_GBOY)
            hc = colorbar;
            set(hc, 'Color', [1 1 1], 'FontName', 'cambria', 'box', 'off','Location', 'southoutside', 'FontWeight', 'bold')   %color bar in white
    
%-----------------------------------------------------------------------------------------------------------------------------            
            
            % Right hemisphere internal side
            cfg              = [];
            cfg.method       = 'surface';
            cfg.funcolorlim  = [min_all abs(min_all)]
            cfg.funparameter = 'pow';
            cfg.camlight     = 'yes';
            cfg.colorbar = 'no'
            cfg.funcolormap   =  'jet';            
            cfg.surffile     = ['C:\Users\Ness\Documents\Toolboxes\fieldtrip-20170116\template\anatomy\surface_white_right.mat']
            ft_sourceplot(cfg, Beta_cond_plot_intp);
            view ([-90 0]);
            set(gcf,'color','black')
            title(['Right hemisphere ' cond],'Color', 'w')
%             set(gca, 'CLim', [2.1 4.8])
%             colormap(mymap_GBOY)
            hc = colorbar;
            set(hc, 'Color', [1 1 1], 'FontName', 'cambria', 'box', 'off','Location', 'southoutside', 'FontWeight', 'bold')   %color bar in white
  
            
            % Right hemisphere external side
            cfg              = [];
            cfg.method       = 'surface';
            cfg.funcolorlim  = [min_all abs(min_all)]
            cfg.funparameter = 'pow';
            cfg.camlight     = 'yes';
            cfg.colorbar = 'no'
            cfg.funcolormap   =  'jet';            
            cfg.surffile     = ['C:\Users\Ness\Documents\Toolboxes\fieldtrip-20170116\template\anatomy\surface_white_right.mat']
            ft_sourceplot(cfg, Beta_cond_plot_intp);
            view ([90 0]);
            set(gcf,'color','black')
            title(['Right hemisphere ' cond],'Color', 'w')
%             set(gca, 'CLim', [2.1 4.8])
%             colormap(mymap_GBOY)
            hc = colorbar;
            set(hc, 'Color', [1 1 1], 'FontName', 'cambria', 'box', 'off','Location', 'southoutside', 'FontWeight', 'bold')   %color bar in white
 
 %-----------------------------------------------------------------------------------------------------------------------------            
 
            % Both hemisphere 
            cfg              = [];
            cfg.method       = 'surface';
            cfg.funcolorlim  = [min_all abs(min_all)]
            cfg.funparameter = 'pow';
            cfg.camlight     = 'yes';
            cfg.funcolormap   =  'jet';    
%             cfg.funcolorlim  = [0 5]
%             cfg.opacitylim  = [0 2]
%             cfg.opacitymap  = 'rampup'
            cfg.surffile     = ['C:\Users\Ness\Documents\Toolboxes\fieldtrip-20170116\template\anatomy\surface_white_both.mat']
            ft_sourceplot(cfg, Beta_cond_plot_intp);
            view ([0 90]);
            set(gcf,'color','black')
            title(['Both hemisphere ' cond],'Color', 'w')
%             set(gca, 'CLim', [2.1 4.8])
%             colormap(mymap_GBOY)
            hc = colorbar;
            set(hc, 'Color', [1 1 1], 'FontName', 'cambria', 'box', 'off','Location', 'southoutside', 'FontWeight', 'bold')   %color bar in white
 
            
