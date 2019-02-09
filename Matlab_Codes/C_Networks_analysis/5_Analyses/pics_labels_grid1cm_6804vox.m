%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% To visualise the results of the GLM
%%%%%
%%%%% -> Display the betas with the highest value (=pics) 
%%%%%    of the componants of interest 
%%%%%    on a standard brain MRI
%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialise

path2ft = '/lena13/home_users/users/hazem/fieldtrip-20161231' 
%path2ft = 'C:\Users\Ness\Documents\Toolboxes\fieldtrip-20170116'

% read the atlas
atlas = ft_read_atlas([ path2ft filesep 'template' filesep 'atlas' filesep 'aal' filesep 'ROI_MNI_V4.nii']) %, filename of atlas to use, see FT_READ_ATLAS

% load the template sourcemodel with the resolution you need (i.e. the resolution you used in your beamformer grid)
load([ path2ft filesep 'template' filesep 'sourcemodel' filesep 'standard_sourcemodel3d10mm.mat'])

%% Create index of voxels for each region of the aal

% Interpolate the atlas and the grid source model of 1cm with ft_sourceinterpolate: 
cfg = []; 
cfg.interpmethod = 'nearest'; 
cfg.parameter = 'tissue'; 
sourcemodel2 = ft_sourceinterpolate(cfg, atlas, sourcemodel); 

% Index of voxels corresponding to each of the 116 labels of the AAL
nbVox = 0
for lab_AAL = 1:116
indx{lab_AAL, 1} = sourcemodel2.tissuelabel(lab_AAL);    
indx{lab_AAL, 2} = find(sourcemodel2.tissue==lab_AAL);
nbVox = nbVox + length(indx{lab_AAL,2});  %brain volum in cm3 occupied by AAL = 1459
end

%% Compute Mean Beta in each cond and each aal ROI for each suj 

%Compute Beta for each ind and each cond , for each ROI (mean of ROI voxels of 1cm3)
path_file = '/pclxserver/raid11/data/SELFI/Networks_Group_comparison/alpha/ICAdecomp/ncomps_15/icasso_nica15_corr25/'
load('/pclxserver/raid11/data/SELFI/Networks_Group_comparison/alpha/SourceData/alpha_betas_cond1_cond2.mat')
%load('E:\Autre\to_send\alpha_betas_cond1_cond2.mat')

%Betas of the comp of interest
beta_comp15.cond1 = squeeze(betas.cond1(12, :, :)); % comp15
beta_comp15.cond2 = squeeze(betas.cond2(12, :, :));
clear betas

%Mean Beta for each cond/ROI/suj
for i = 1:length(indx(:, 1)) % pour tous les ROIs  
indx_betas = indx{i, 2}(:); %index des voxels du ROI

for n = 1: 19 % nb suj
indx{i, 3}(n, 1) = {mean(beta_comp15.cond1(indx_betas, n))}; % beta_cond1-ROIAAL
indx{i, 4}(n, 1) = {mean(beta_comp15.cond2(indx_betas, n))}; % beta_cond2-ROIAAL
end

end

%Beta for each cond/suj in a specific voxel (here : 5212)
indx_betas = 5212; %index des voxels du ROI

for n = 1: 19 % nb suj
indx_max{n, 1} = indx_betas
indx_max{n, 2} = beta_comp15.cond1(indx_betas, n); % beta_cond1
indx_max{n, 3} = beta_comp15.cond2(indx_betas, n); % beta_cond2
end

%% Create index of voxel in the cluster (t>2.4450) and cluster pics (t>3.1966)

%files with the cluster pvalues
load([path_file 'stat_comp_5000rands_cleancorr_clusteralpha025.mat'])
%load('E:\data_SELFI\last_analyses\alphaicasso_nica15_corr25\stat_comp_5000rands_cleancorr_clusteralpha025.mat')
tval = stat_comp.comp11.stat;

%List of voxels (on the 6804) included in the cluster, and of the pics of the cluster 
indx_cluster{:, 1} = find((abs(tval(:,1)) > 2.4450) & (~isnan(tval(:,1)))); %1122
indx_pics{:, 1} = find((abs(tval(:,1)) > 3.1966) & (~isnan(tval(:,1)))); %231

%List of voxels not NAN 
indx_brain{:, 1} = find((~isnan(tval(:,1)))); %4933

%%  For each aal ROI: Create index of voxel of the ROI in the cluster (t>2.4450) and cluster pics (t>3.1966)

%Mean Beta for each cond/ROI/suj
Reshape_sourcemodel2.tissue = reshape(sourcemodel2.tissue, [18*21*18, 1]);

for lab_AAL = 1:116 % pour tous les ROIs         
size_ROI = length(indx{lab_AAL, 2}(:)); %taille du ROI pour verif
indx_betasclusterROI = find((abs(tval(:,1)) > 2.4450) & (~isnan(tval(:,1))) & (Reshape_sourcemodel2.tissue==lab_AAL));
indx_betaspicsROI = find((abs(tval(:,1)) > 3.1966) & (~isnan(tval(:,1))) & (Reshape_sourcemodel2.tissue==lab_AAL));

for n = 1: 19 % nb suj
indx{lab_AAL, 5}(n, 1) = {mean(beta_comp15.cond1(indx_betasclusterROI, n))}; % beta_cond1 - ROI & (t>2.4450)
indx{lab_AAL, 6}(n, 1) = {mean(beta_comp15.cond2(indx_betasclusterROI, n))}; % beta_cond2 - ROI & (t>2.4450)
indx{lab_AAL, 8}(n, 1) = {mean(beta_comp15.cond1(indx_betaspicsROI, n))}; % beta_cond1 - ROI & (t>3.1966)
indx{lab_AAL, 9}(n, 1) = {mean(beta_comp15.cond2(indx_betaspicsROI, n))}; % beta_cond2 - ROI & (t>3.1966)
end
indx{lab_AAL, 7} = length(indx_betasclusterROI);  % nbr de points to compute mean beta- ROI & (t>2.4450)
indx{lab_AAL, 10} = length(indx_betaspicsROI); % nbr de points to compute mean beta- ROI & (t>3.1966)
end

%% 117: Mean for 'cing ant', 'cing med' and 'front sup med' right taken as one region
indx{117, 1} = 'Cing_Ant-Cing_Med-Front_Sup_Medial-R'
indx{117, 2} = find((Reshape_sourcemodel2.tissue==32) | (Reshape_sourcemodel2.tissue==34)| (Reshape_sourcemodel2.tissue==24));

size_117 = (length(indx{32, 2}(:)) + length(indx{34, 2}(:)) + length(indx{24, 2}(:))); %for verif size index

indx_betas = indx{117, 2}(:);
indx_betasclusterROI = find((abs(tval(:,1)) > 2.4450) & (~isnan(tval(:,1))) & ((Reshape_sourcemodel2.tissue==32) | (Reshape_sourcemodel2.tissue==34)| (Reshape_sourcemodel2.tissue==24)));
indx_betaspicsROI = find((abs(tval(:,1)) > 3.1966) & (~isnan(tval(:,1))) & ((Reshape_sourcemodel2.tissue==32) | (Reshape_sourcemodel2.tissue==34)| (Reshape_sourcemodel2.tissue==24)));

for n = 1: 19 % nb suj
indx{117, 3}(n, 1) = {mean(beta_comp15.cond1(indx_betas, n))}; % beta_cond1-ROIAAL
indx{117, 4}(n, 1) = {mean(beta_comp15.cond2(indx_betas, n))} ;% beta_cond2-ROIAAL
indx{117, 5}(n, 1) = {mean(beta_comp15.cond1(indx_betasclusterROI, n))}; % beta_cond1 - ROI & (t>2.4450)
indx{117, 6}(n, 1) = {mean(beta_comp15.cond2(indx_betasclusterROI, n))} ;% beta_cond2 - ROI & (t>2.4450)
indx{117, 8}(n, 1) = {mean(beta_comp15.cond1(indx_betaspicsROI, n))} ;% beta_cond1 - ROI & (t>3.1966)
indx{117, 9}(n, 1) = {mean(beta_comp15.cond2(indx_betaspicsROI, n))}; % beta_cond2 - ROI & (t>3.1966)
end
indx{117, 7} = length(indx_betasclusterROI);  % nbr de points to compute mean beta- ROI & (t>2.4450)
indx{117, 10} = length(indx_betaspicsROI) ;% nbr de points to compute mean beta- ROI & (t>3.1966)

%% Mean for several regions taken as one region !!!to adapt as a fonction of the number of regions
nb_line = length(indx) +1
name_regions = 'Frontal R : sup, sup orb, mid, mid orb, inf oper, inf tri, inf orb'
aal_number_regions = [4 6 8 10 12 14 16]

indx{nb_line, 1} = name_regions
indx{nb_line, 2} = find((Reshape_sourcemodel2.tissue==aal_number_regions(1)) | (Reshape_sourcemodel2.tissue==aal_number_regions(2)) | (Reshape_sourcemodel2.tissue==aal_number_regions(3)) | (Reshape_sourcemodel2.tissue==aal_number_regions(4)) | (Reshape_sourcemodel2.tissue==aal_number_regions(5)) | (Reshape_sourcemodel2.tissue==aal_number_regions(6)) | (Reshape_sourcemodel2.tissue==aal_number_regions(7)));
%size_region = (length(indx{aal_number_regions(1), 2}(:)) + length(indx{aal_number_regions(2), 2}(:)) + length(indx{aal_number_regions(3), 2}(:))); %for verif size index

indx_betas = indx{nb_line, 2}(:);
indx_betasclusterROI = find((abs(tval(:,1)) > 2.4450) & (~isnan(tval(:,1))) & ((Reshape_sourcemodel2.tissue==aal_number_regions(1)) | (Reshape_sourcemodel2.tissue==aal_number_regions(2))| (Reshape_sourcemodel2.tissue==aal_number_regions(3)) | (Reshape_sourcemodel2.tissue==aal_number_regions(4)) | (Reshape_sourcemodel2.tissue==aal_number_regions(5)) | (Reshape_sourcemodel2.tissue==aal_number_regions(6)) | (Reshape_sourcemodel2.tissue==aal_number_regions(7))));
indx_betaspicsROI = find((abs(tval(:,1)) > 3.1966) & (~isnan(tval(:,1))) & ((Reshape_sourcemodel2.tissue==aal_number_regions(1)) | (Reshape_sourcemodel2.tissue==aal_number_regions(2))| (Reshape_sourcemodel2.tissue==aal_number_regions(3)) | (Reshape_sourcemodel2.tissue==aal_number_regions(4)) | (Reshape_sourcemodel2.tissue==aal_number_regions(5)) | (Reshape_sourcemodel2.tissue==aal_number_regions(6)) | (Reshape_sourcemodel2.tissue==aal_number_regions(7))));

for n = 1: 19 % nb suj
indx{nb_line, 3}(n, 1) = {mean(beta_comp15.cond1(indx_betas, n))}; % beta_cond1-ROIAAL
indx{nb_line, 4}(n, 1) = {mean(beta_comp15.cond2(indx_betas, n))} ;% beta_cond2-ROIAAL
indx{nb_line, 5}(n, 1) = {mean(beta_comp15.cond1(indx_betasclusterROI, n))}; % beta_cond1 - ROI & (t>2.4450)
indx{nb_line, 6}(n, 1) = {mean(beta_comp15.cond2(indx_betasclusterROI, n))} ;% beta_cond2 - ROI & (t>2.4450)
indx{nb_line, 8}(n, 1) = {mean(beta_comp15.cond1(indx_betaspicsROI, n))} ;% beta_cond1 - ROI & (t>3.1966)
indx{nb_line, 9}(n, 1) = {mean(beta_comp15.cond2(indx_betaspicsROI, n))}; % beta_cond2 - ROI & (t>3.1966)
end
indx{nb_line, 7} = length(indx_betasclusterROI);  % nbr de points to compute mean beta- ROI & (t>2.4450)
indx{nb_line, 10} = length(indx_betaspicsROI) ;% nbr de points to compute mean beta- ROI & (t>3.1966)

%% Plot
% 
% template_grid = load(['C:\Users\Ness\Documents\Toolboxes\fieldtrip-20170116\template\sourcemodel\standard_sourcemodel3d10mm.mat']);
% templatedir  = ['C:\Users\Ness\Documents\Toolboxes\fieldtrip-20170116\external\spm8\templates'];
% template = ft_read_mri(fullfile(templatedir, 'T1.nii'));
% template_resliced = ft_volumereslice([], template);  
% 
% % Create a fieldtrip-like data structure with the pvalues
% pval_plot = [];
% pval_plot.dim = template_grid.sourcemodel.dim;
% pval_plot.inside = template_grid.sourcemodel.inside;
% pval_plot.pos = template_grid.sourcemodel.pos;
% pval_plot.method = 'average';
% pval_plot.time{1} = pval_plot2.time %data_hilbert_all.time;
% pval_plot.pow = tval; % here are the p-values for each grid point
% 
% % Mask before interpolation
% pval_plot.pow(:, 1) = pval_plot.pow(:, 1)>2.4450  
% pval_plot.pow(pval_plot.pow>2.4450) == 1
% 
% % Interpolation: Project on standard MRI brain
% cfg             = [];
% cfg.parameter   = 'pow';
% pval_plot_intp  = ft_sourceinterpolate(cfg, pval_plot, template_resliced);
% 
% 
% %Visualisation
%             cfg              = [];
%             cfg.method       = 'surface';
%             cfg.funcolorlim  = 'zeromax' %'minzero'%'maxabs'
%             cfg.funparameter = 'pow';
%             cfg.camlight     = 'yes';
% %             cfg.colorbar = 'no'
%             cfg.funcolormap   =  'hot';             
%             cfg.surffile     = ['C:\Users\Ness\Documents\Toolboxes\fieldtrip-20170116\template\anatomy\surface_white_right.mat']
%             ft_sourceplot(cfg, pval_plot_intp);
%             view ([-90 0]);
%             set(gcf,'color','black')
%             title('Right hemisphere ','Color', 'w')




