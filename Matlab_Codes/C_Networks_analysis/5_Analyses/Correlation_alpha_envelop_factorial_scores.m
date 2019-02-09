%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SCRIPT_Correlation_alpha_envelop_factorial_scores
%%%
%%% For each voxel, and each dimension 
%%% compute the correlation between the mean alpha enveloppe (over the 10min)
%%% and the factorial scores over the participants 
%%% in the high contact condition (19 points), the low contact condition (19 points)
%%% and both condition taken together (38 points)
%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load concatenated source data
freq_band = 'alpha'
Output = ['/pclxserver/raid11/data/SELFI/Networks_Group_comparison/' freq_band '/SourceData'];
load([Output filesep freq_band '_hilbert_across_subjects_across_groups_smooth.mat']);

% Create vector indicating for each data point to which subject it corresponds to
vect_subj = ones(1, length(data_hilbert_all.condition));
ind_subj = 1;
for i=2:length(vect_subj)
    if data_hilbert_all.condition(i) == 11 && data_hilbert_all.condition(i-1)~=11
        ind_subj = ind_subj+1;
    end
    vect_subj(i) = ind_subj;
end

% Create consistent label field
for i=1:size(data_hilbert_all.trial{1}, 1)
    data_hilbert_all.label{i} = ['v' num2str(i)];
end

% % Remove voxels = 0 (outside the brain)
% nnz_vox = find(data_hilbert_all.trial{1}(:, 1) ~= 0);
% cfg = [];
% cfg.channel = nnz_vox;
% data_hilbert_all_nnz = ft_selectdata(cfg, data_hilbert_all);
% data_hilbert_all_nnz.pos = data_hilbert_all_nnz.pos(nnz_vox, :);

%% Compute alpha mean value in each voxel, for each sub, in each cond (38 points)

vect_condH = (data_hilbert_all.condition == 11) | (data_hilbert_all.condition == 12);
vect_condL = (data_hilbert_all.condition == 21) | (data_hilbert_all.condition == 22);

    for ngrid = 1:size(data_hilbert_all.trial{1},1)
        for subj = 1:19
Vect_condH_subj = (vect_condH & vect_subj == subj) ;         
Vect_condL_subj = (vect_condL & vect_subj == subj) ;       

alphamean.value(ngrid, subj) = mean(data_hilbert_all.trial{1}(ngrid, Vect_condH_subj));
alphamean.subj(ngrid,subj) = subj;

alphamean.value(ngrid, (subj+19)) = mean(data_hilbert_all.trial{1}(ngrid, Vect_condL_subj));
alphamean.subj(ngrid, (subj+19)) = subj;
        end
    end
    
  alphamean.cond(ngrid, 1:19) = 1;
  alphamean.cond(ngrid, 20:38) = 2;

 save([Output filesep 'alphamean'], 'alphamean', '-v7.3')
 
  %% Create a matrix with the factorial scores (for each dim: s1H, s2H, s3H,...., s1L, s2L....)
  
  factorialscores.MFA1.dim1 = zeros(1 , 38)
  factorialscores.MFA1.dim2 = zeros(1 , 38)
  factorialscores.MFA1.dim3 = zeros(1 , 38)
  
  factorialscores.MFA3.dim3 = zeros(1 , 38)
  
  %Paste from excel the factorial scores
  
  factorialscores.MFA1.dim1 = factorialscores.MFA1.dim1(1, 1:38)
  factorialscores.MFA1.dim2 = factorialscores.MFA1.dim2(1, 1:38)
  factorialscores.MFA1.dim3 = factorialscores.MFA1.dim3(1, 1:38)
  
  factorialscores.MFA3.dim3 = factorialscores.MFA3.dim3(1, 1:38)
  
  save([Output filesep 'factorialscores_MFA1_MFA3'], 'factorialscores', '-v7.3')
  
  %% Compute the correlation between the mean alpha value and the factorial scores
  
   for ngrid = 1:size(data_hilbert_all.trial{1},1)
    for dim = 1:3
  [R, P] = corrcoef(alphamean.value(ngrid, 1:19), factorialscores.MFA1.(['dim' num2str(dim)])(1, 1:19))  
  corr.(['dim' num2str(dim)]).High.R(ngrid, 1) = R(1, 2)
  corr.(['dim' num2str(dim)]).High.P(ngrid, 1) = P(1, 2)
  
  [R, P] = corrcoef(alphamean.value(ngrid, 20:38), factorialscores.MFA1.(['dim' num2str(dim)])(1, 20:38))  
  corr.(['dim' num2str(dim)]).Low.R(ngrid, 1) = R(1, 2)
  corr.(['dim' num2str(dim)]).Low.P(ngrid, 1) = P(1, 2)
  
  [R, P] = corrcoef(alphamean.value(ngrid, 1:38), factorialscores.MFA1.(['dim' num2str(dim)])(1, 1:38))  
  corr.(['dim' num2str(dim)]).HighLow.R(ngrid, 1) = R(1, 2)
  corr.(['dim' num2str(dim)]).HighLow.P(ngrid, 1) = P(1, 2)
  
    end
    end
  
    save([Output filesep 'corr_MFA1'], 'corr', '-v7.3')
    
    for ngrid = 1:size(data_hilbert_all.trial{1},1)
  dim = 3
  
  [R, P] = corrcoef(alphamean.value(ngrid, 1:19), factorialscores.MFA3.(['dim' num2str(dim)])(1, 1:19))  
  corr.(['dim' num2str(dim)]).High.R(ngrid, 1) = R(1, 2)
  corr.(['dim' num2str(dim)]).High.P(ngrid, 1) = P(1, 2)
  
  [R, P] = corrcoef(alphamean.value(ngrid, 20:38), factorialscores.MFA3.(['dim' num2str(dim)])(1, 20:38))  
  corr.(['dim' num2str(dim)]).Low.R(ngrid, 1) = R(1, 2)
  corr.(['dim' num2str(dim)]).Low.P(ngrid, 1) = P(1, 2)
  
  [R, P] = corrcoef(alphamean.value(ngrid, 1:38), factorialscores.MFA3.(['dim' num2str(dim)])(1, 1:38))  
  corr.(['dim' num2str(dim)]).HighLow.R(ngrid, 1) = R(1, 2)
  corr.(['dim' num2str(dim)]).HighLow.P(ngrid, 1) = P(1, 2)
  
    end
  
    save([Output filesep 'corr_MFA3'], 'corr', '-v7.3')
  
   %% For each dimension: Plot the correlation value of the voxel with p<0.05
   
path2ft = '/lena13/home_users/users/hazem/fieldtrip-20161231'   
   
template_grid = load([path2ft filesep 'template' filesep 'sourcemodel' filesep 'standard_sourcemodel3d10mm.mat']);
templatedir  = [path2ft filesep 'external' filesep 'spm8' filesep 'templates'];
template = ft_read_mri(fullfile(templatedir, 'T1.nii'));
template_resliced = ft_volumereslice([], template);  
      
%load([Output filesep 'corr_MFA1'])

load([Output filesep 'corr_MFA3'])

    for dim = 1:3
       
%cond = 'High'
%cond = 'Low'
%cond = 'HighLow'

% Create a fieldtrip-like data structure with the correlation values
correlation_plot = [];
correlation_plot.dim = template_grid.sourcemodel.dim;
correlation_plot.inside = template_grid.sourcemodel.inside;
correlation_plot.pos = template_grid.sourcemodel.pos;
correlation_plot.method = 'average';
correlation_plot.time = data_hilbert_all.time;
correlation_plot.avg.pow = corr.(['dim' num2str(dim)]).(cond).R; % here are the correltaion values for each grid point

% Mask: corr = 0 for non significant P values
 correlation_plot.avg.pow(corr.(['dim' num2str(dim)]).(cond).P > 0.05) = 0         

% Interpolation: Project computed sources on standard MRI brain
cfg            = [];
cfg.parameter = 'avg.pow';
correlation_plot_intp  = ft_sourceinterpolate(cfg, correlation_plot, template_resliced);


            % Left hemisphere internal side
            cfg              = [];
            cfg.method       = 'surface';
            cfg.funcolorlim  = 'maxabs' %'maxabs' %'minzero'%'maxabs'
            cfg.funparameter = 'pow';
            cfg.camlight     = 'yes';
             cfg.colorbar = 'no'
            cfg.funcolormap   =  'jet';            
            cfg.surffile     = [path2ft filesep 'template' filesep 'anatomy' filesep 'surface_white_left.mat']
%               correlation_plot_intp.mask(:, 1) = ~(correlation_plot_intp.pow(:, 1) == 0) %Mask not significant correlations
%               cfg.maskparameter = 'mask' 
            ft_sourceplot(cfg, correlation_plot_intp);
            view ([90 0]);
            set(gcf,'color','black')
            title(['Left hemisphere ' cond ' dim ' num2str(dim)],'Color', 'w')
            hc = colorbar;
            set(hc, 'Color', [1 1 1], 'FontName', 'cambria', 'box', 'off','Location', 'southoutside', 'FontWeight', 'bold')   %color bar in white
                         
            % Left hemisphere external side
            cfg              = [];
            cfg.method       = 'surface';
            cfg.funcolorlim  = 'maxabs' 
            cfg.funparameter = 'pow';
            cfg.camlight     = 'yes';
             cfg.colorbar = 'no'
            cfg.funcolormap   =  'jet';          
            cfg.surffile     = [path2ft filesep 'template' filesep 'anatomy' filesep 'surface_white_left.mat']
            ft_sourceplot(cfg, correlation_plot_intp);
            view ([-90 0]);
            set(gcf,'color','black')
            title(['Left hemisphere ' cond ' dim ' num2str(dim)],'Color', 'w')
            hc = colorbar;
            set(hc, 'Color', [1 1 1], 'FontName', 'cambria', 'box', 'off','Location', 'southoutside', 'FontWeight', 'bold')   %color bar in white
    
%-----------------------------------------------------------------------------------------------------------------------------            
            
            % Right hemisphere internal side
            cfg              = [];
            cfg.method       = 'surface';
            cfg.funcolorlim  = 'maxabs'
            cfg.funparameter = 'pow';
            cfg.camlight     = 'yes';
            cfg.colorbar = 'no'
            cfg.funcolormap   =  'jet';            
            cfg.surffile     = [path2ft filesep 'template' filesep 'anatomy' filesep 'surface_white_right.mat']
            ft_sourceplot(cfg, correlation_plot_intp);
            view ([-90 0]);
            set(gcf,'color','black')
            title(['Right hemisphere ' cond ' dim ' num2str(dim)],'Color', 'w')
            hc = colorbar;
            set(hc, 'Color', [1 1 1], 'FontName', 'cambria', 'box', 'off','Location', 'southoutside', 'FontWeight', 'bold')   %color bar in white
  
            
            % Right hemisphere external side
            cfg              = [];
            cfg.method       = 'surface';
            cfg.funcolorlim  = 'maxabs'
            cfg.funparameter = 'pow';
            cfg.camlight     = 'yes';
             cfg.colorbar = 'no'
            cfg.funcolormap   =  'jet';            
            cfg.surffile     = [path2ft filesep 'template' filesep 'anatomy' filesep 'surface_white_right.mat']
            ft_sourceplot(cfg, correlation_plot_intp);
            view ([90 0]);
            set(gcf,'color','black')
            title(['Right hemisphere ' cond ' dim ' num2str(dim)],'Color', 'w')
            hc = colorbar;
            set(hc, 'Color', [1 1 1], 'FontName', 'cambria', 'box', 'off','Location', 'southoutside', 'FontWeight', 'bold')   %color bar in white
 
 %-----------------------------------------------------------------------------------------------------------------------------            
 
            % Both hemisphere 
            cfg              = [];
            cfg.method       = 'surface';
            cfg.funcolorlim  = 'maxabs'
            cfg.funparameter = 'pow';
            cfg.camlight     = 'yes';
            cfg.funcolormap   =  'jet';  
             cfg.colorbar = 'no'
%             cfg.funcolorlim  = [0 5]
%             cfg.opacitylim  = [0 2]
%             cfg.opacitymap  = 'rampup'
            cfg.surffile     = [path2ft filesep 'template' filesep 'anatomy' filesep 'surface_white_both.mat']
            ft_sourceplot(cfg, correlation_plot_intp);
            view ([0 90]);
            set(gcf,'color','black')
            title(['Both hemisphere ' cond ' dim ' num2str(dim)],'Color', 'w')
            hc = colorbar;
            set(hc, 'Color', [1 1 1], 'FontName', 'cambria', 'box', 'off','Location', 'southoutside', 'FontWeight', 'bold')   %color bar in white
 
            
    end
  
  for n = 1:length(Table(:, 1))
  Table_indvalues_significant.alphamean(n, 1:38) = alphamean.value(Table(n, 1), 1:38)
  
  end
  
  Table_indvalues_significant.factorialscores_Mfa3_dim3(1, 1:38) = factorialscores.MFA3.dim3(1, 1:38)
  
  Table_indvalues_significant.list_corrval = Table
  
