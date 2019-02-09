%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% To visualise the results of the GLM
%%%%%
%%%%% -> Display the betas
%%%%%    of the componants of interest 
%%%%%    on a surface view of the brain
%%%%%    using own made color map with more marked contrast
%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-----------------------------------------------------------------------------
% load('/pclxserver/raid11/data/SELFI/archive_sourceplot_surface/pval_plot_int.mat')
%-----------------------------------------------------------------------------

load('/pclxserver/raid11/data/SELFI//Networks_Group_comparison/alpha/SourceData/alpha_hilbert_across_subjects_across_groups_smooth.mat');
load('/pclxserver/raid11/data/SELFI/Networks_Group_comparison/alpha/ICAdecomp/ncomps_15/icasso_nica15_corr25/stat_comp_5000rands_cleancorr_clusteralpha025.mat')
plot_tIC = stat_comp.comp11.stat;

% Option with MRI from one subject (clearer images and allows to apply mask on image)
% template_grid = load('/pclxserver/raid11/data/SELFI/selfi01_s01/Anatomy/grid1cm.mat');
% template_grid.sourcemodel = template_grid.grid1cm;
% template_grid = rmfield(template_grid, 'grid1cm');
% template_resliced = load('/pclxserver/raid11/data/SELFI/selfi01_s01/Anatomy/mri_realigned_resliced_meters.mat');
% template_resliced = template_resliced.mri_realigned_resliced_meters;
% load('/pclxserver/raid11/data/SELFI/selfi01_s01/Anatomy/segmentedmri.mat');
% template_resliced.anatomy = template_resliced.anatomy .* segmentedmri.brain;

template_grid = load(['C:\Users\Ness\Documents\Toolboxes\fieldtrip-20170116\template\sourcemodel\standard_sourcemodel3d10mm.mat']);
templatedir  = ['C:\Users\Ness\Documents\Toolboxes\fieldtrip-20170116\external\spm8\templates'];
template = ft_read_mri(fullfile(templatedir, 'T1.nii'));
template_resliced = ft_volumereslice([], template);  

% template_grid = load(['/lena13/home_users/users/hazem/fieldtrip-20161231/template/sourcemodel/standard_sourcemodel3d10mm.mat']);
% templatedir  = ['/lena13/home_users/users/hazem/fieldtrip-20161231/external/spm8/templates'];
% template = ft_read_mri(fullfile(templatedir, 'T1.nii'));
% template_resliced = ft_volumereslice([], template);  

 % Create a fieldtrip-like data structure with the pvalues
pval_plot = [];
pval_plot.dim = template_grid.sourcemodel.dim;
pval_plot.inside = template_grid.sourcemodel.inside;
pval_plot.pos = template_grid.sourcemodel.pos;
pval_plot.method = 'average';
pval_plot.time{1} = pval_plot2.time %data_hilbert_all.time;
pval_plot.avg.pow = plot_tIC; % here are the p-values for each grid point
            
% Interpolation: Project on standard MRI brain
cfg             = [];
cfg.parameter   = 'pow';
pval_plot_intp  = ft_sourceinterpolate(cfg, pval_plot, template_resliced);
            
%---------------------------------------------------------------------------
% tval_0.05= 2.1009
% tval_0.025= 2.4450
% tval_0.01= 2.8784
% tval_0.005= 3.1966

mymap = mymap_5hotColor
mymap.color = mymap_red
            % Left hemisphere internal side
            cfg              = [];
            cfg.method       = 'surface';
            cfg.funcolorlim  = 'zeromax' %'minzero'%'maxabs'
            cfg.funparameter = 'pow';
            cfg.camlight     = 'yes';
            cfg.colorbar = 'no'
            cfg.funcolormap   =  'hot';            
            pval_plot_intp.mask(:, 1) = (pval_plot_intp.pow(:, 1)>2.4450) 
            cfg.maskparameter = 'mask'
            cfg.surffile     = ['C:\Users\Ness\Documents\Toolboxes\fieldtrip-20170116\template\anatomy\surface_white_left.mat']
            ft_sourceplot(cfg, pval_plot_intp);
            view ([90 0]);
            set(gcf,'color','black')
            title('Left hemisphere ','Color', 'w')
            set(gca, 'CLim', [2 4.7805])
            colormap(mymap.color)
            hc = colorbar;
%             set(hc, 'Color', [1 1 1], 'FontName', 'cambria', 'box', 'off','Location', 'southoutside', 'FontWeight', 'bold', 'Limits', [2.1  4.8], 'Ticks', [2.1, 2.9, 3.2, 3.9, 4.8], 'Ticklabels', {'p=0.05','p=0.01','p=0.005', 'p=0.001', 'p=0.00016'})   %color bar in white
            set(hc, 'Color', [1 1 1], 'FontName', 'cambria', 'box', 'off','Location', 'southoutside', 'FontWeight', 'bold', 'Limits', [2.4450  4.7805], 'Ticks', [mymap.Ticks(1, :), mymap.Ticks(2, :), mymap.Ticks(3, :), mymap.Ticks(4, :), mymap.Ticks(5, :), mymap.Ticks(6, :)], 'Ticklabels', {mymap.Index(1, :), mymap.Index(2, :), mymap.Index(3, :), mymap.Index(4, :), mymap.Index(5, :), mymap.Index(6, :)})   %color bar in white
                          
            % Left hemisphere external side
            cfg              = [];
            cfg.method       = 'surface';
            cfg.funcolorlim  = 'zeromax' %'minzero'%'maxabs'
            cfg.funparameter = 'pow';
            cfg.camlight     = 'yes';
            cfg.colorbar = 'no'
            cfg.funcolormap   =  'hot';            
            pval_plot_intp.mask(:, 1) = (pval_plot_intp.pow(:, 1)>2.4450)  
            cfg.maskparameter = 'mask'
            cfg.surffile     = ['C:\Users\Ness\Documents\Toolboxes\fieldtrip-20170116\template\anatomy\surface_white_left.mat']
            ft_sourceplot(cfg, pval_plot_intp);
            view ([-90 0]);
            set(gcf,'color','black')
            title('Left hemisphere ','Color', 'w')
            set(gca, 'CLim', [2 4.7805])
            colormap(mymap.color)
            hc = colorbar;
%             set(hc, 'Color', [1 1 1], 'FontName', 'cambria', 'box', 'off','Location', 'southoutside', 'FontWeight', 'bold', 'Limits', [2.1  4.8], 'Ticks', [2.1, 2.9, 3.2, 3.9, 4.8], 'Ticklabels', {'p=0.05','p=0.01','p=0.005', 'p=0.001', 'p=0.00016'})   %color bar in white
            set(hc, 'Color', [1 1 1], 'FontName', 'cambria', 'box', 'off','Location', 'southoutside', 'FontWeight', 'bold', 'Limits', [2.4450  4.7805], 'Ticks', [mymap.Ticks(1, :), mymap.Ticks(2, :), mymap.Ticks(3, :), mymap.Ticks(4, :), mymap.Ticks(5, :), mymap.Ticks(6, :)], 'Ticklabels', {mymap.Index(1, :), mymap.Index(2, :), mymap.Index(3, :), mymap.Index(4, :), mymap.Index(5, :), mymap.Index(6, :)})   %color bar in white
              
%-----------------------------------------------------------------------------------------------------------------------------            
            
            % Right hemisphere internal side
            cfg              = [];
            cfg.method       = 'surface';
            cfg.funcolorlim  = 'zeromax' %'minzero'%'maxabs'
            cfg.funparameter = 'pow';
            cfg.camlight     = 'yes';
            cfg.colorbar = 'no'
            cfg.funcolormap   =  'hot';            
            pval_plot_intp.mask(:, 1) = (pval_plot_intp.pow(:, 1)>2.4450) 
            cfg.maskparameter = 'mask'
            cfg.surffile     = ['C:\Users\Ness\Documents\Toolboxes\fieldtrip-20170116\template\anatomy\surface_white_right.mat']
            ft_sourceplot(cfg, pval_plot_intp);
            view ([-90 0]);
            set(gcf,'color','black')
            title('Right hemisphere ','Color', 'w')
            set(gca, 'CLim', [2 4.7805])
            colormap(mymap.color)
            hc = colorbar;
%             set(hc, 'Color', [1 1 1], 'FontName', 'cambria', 'box', 'off','Location', 'southoutside', 'FontWeight', 'bold', 'Limits', [2.1  4.8], 'Ticks', [2.1, 2.9, 3.2, 3.9, 4.8], 'Ticklabels', {'p=0.05','p=0.01','p=0.005', 'p=0.001', 'p=0.00016'})   %color bar in white
            set(hc, 'Color', [1 1 1], 'FontName', 'cambria', 'box', 'off','Location', 'southoutside', 'FontWeight', 'bold', 'Limits', [2.4450  4.7805], 'Ticks', [mymap.Ticks(1, :), mymap.Ticks(2, :), mymap.Ticks(3, :), mymap.Ticks(4, :), mymap.Ticks(5, :), mymap.Ticks(6, :)], 'Ticklabels', {mymap.Index(1, :), mymap.Index(2, :), mymap.Index(3, :), mymap.Index(4, :), mymap.Index(5, :), mymap.Index(6, :)})   %color bar in white
              
            
            % Right hemisphere external side
            cfg              = [];
            cfg.method       = 'surface';
            cfg.funcolorlim  = 'zeromax' %'minzero'%'maxabs'
            cfg.funparameter = 'pow';
            cfg.camlight     = 'yes';
            cfg.colorbar = 'no'
            cfg.funcolormap   =  'hot';            
            pval_plot_intp.mask(:, 1) = (pval_plot_intp.pow(:, 1)>2.4450) 
            cfg.maskparameter = 'mask'
            cfg.surffile     = ['C:\Users\Ness\Documents\Toolboxes\fieldtrip-20170116\template\anatomy\surface_white_right.mat']
            ft_sourceplot(cfg, pval_plot_intp);
            view ([90 0]);
            set(gcf,'color','black')
            title('Right hemisphere ','Color', 'w')
            set(gca, 'CLim', [2 4.7805])
            colormap(mymap.color)
            hc = colorbar;
%             set(hc, 'Color', [1 1 1], 'FontName', 'cambria', 'box', 'off','Location', 'southoutside', 'FontWeight', 'bold', 'Limits', [2.1  4.8], 'Ticks', [2.1, 2.9, 3.2, 3.9, 4.8], 'Ticklabels', {'p=0.05','p=0.01','p=0.005', 'p=0.001', 'p=0.00016'})   %color bar in white
            set(hc, 'Color', [1 1 1], 'FontName', 'cambria', 'box', 'off','Location', 'southoutside', 'FontWeight', 'bold', 'Limits', [2.4450  4.7805], 'Ticks', [mymap.Ticks(1, :), mymap.Ticks(2, :), mymap.Ticks(3, :), mymap.Ticks(4, :), mymap.Ticks(5, :), mymap.Ticks(6, :)], 'Ticklabels', {mymap.Index(1, :), mymap.Index(2, :), mymap.Index(3, :), mymap.Index(4, :), mymap.Index(5, :), mymap.Index(6, :)})   %color bar in white
              
            
 %-----------------------------------------------------------------------------------------------------------------------------            
 
            % Both hemisphere 
            cfg              = [];
            cfg.method       = 'surface';
            cfg.funcolorlim  = 'zeromax' %'minzero'%'maxabs'
            cfg.funparameter = 'pow';
            cfg.camlight     = 'yes';
            cfg.funcolormap   =  'hot';    
            pval_plot_intp.mask(:, 1) = (pval_plot_intp.pow(:, 1)>2.4450) 
            cfg.maskparameter = 'mask'
%             cfg.funcolorlim  = [0 5]
%             cfg.opacitylim  = [0 2]
%             cfg.opacitymap  = 'rampup'
            cfg.surffile     = ['C:\Users\Ness\Documents\Toolboxes\fieldtrip-20170116\template\anatomy\surface_white_both.mat']
            ft_sourceplot(cfg, pval_plot_intp);
            view ([0 90]);
            set(gcf,'color','black')
            title(['Both hemisphere'] ,'Color', 'w')
            set(gca, 'CLim', [2 4.7805])
            colormap(mymap.color)
            hc = colorbar;
%             set(hc, 'Color', [1 1 1], 'FontName', 'cambria', 'box', 'off','Location', 'southoutside', 'FontWeight', 'bold', 'Limits', [2.1  4.8], 'Ticks', [2.1, 2.9, 3.2, 3.9, 4.8], 'Ticklabels', {'p=0.05','p=0.01','p=0.005', 'p=0.001', 'p=0.00016'})   %color bar in white
            set(hc, 'Color', [1 1 1], 'FontName', 'cambria', 'box', 'off','Location', 'southoutside', 'FontWeight', 'bold', 'Limits', [2.4450  4.7805], 'Ticks', [mymap.Ticks(1, :), mymap.Ticks(2, :), mymap.Ticks(3, :), mymap.Ticks(4, :), mymap.Ticks(5, :), mymap.Ticks(6, :)], 'Ticklabels', {mymap.Index(1, :), mymap.Index(2, :), mymap.Index(3, :), mymap.Index(4, :), mymap.Index(5, :), mymap.Index(6, :)})   %color bar in white
              
            
            