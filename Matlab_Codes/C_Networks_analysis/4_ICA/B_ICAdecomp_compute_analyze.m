%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SCRIPT_ICAdecomp_compute_analyze
%%%
%%% Compute an ICA on the concatenated data accross subjects
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc

%% Init
freq_band = 'theta_b'%'alpha'; %'beta' %'theta_b'
compute_ICA = 'no';  % 'yes': compute ICA, 'no': load existing ICA
numcomponent = 15;
threshold_corr = 0.25;

%% Set path
%  Adapt the path to fieldtrip if user ~= Mariana  :-)
Folder_path = ['/pclxserver/raid11/data/SELFI/Networks_Group_comparison/' freq_band];
%path2ft = '~/../datalinks/INTRA_EMO/INTRA_EMO_Mariana/Scripts/Toolboxes/fieldtrip-20180129';% 
path2ft = '/lena13/home_users/users/hazem/fieldtrip-20161231';
addpath(path2ft);
ft_defaults

%% Load concatenated source data
disp('Loading data')
load([Folder_path  '/SourceData/' freq_band '_hilbert_across_subjects_across_groups_smooth.mat']);

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

% Remove voxels = 0 (outside the brain)
nnz_vox = find(data_hilbert_all.trial{1}(:, 1) ~= 0);
cfg = [];
cfg.channel = nnz_vox;
data_hilbert_all_nnz = ft_selectdata(cfg, data_hilbert_all);
data_hilbert_all_nnz.pos = data_hilbert_all_nnz.pos(nnz_vox, :);

%% Compute ICA with fastica via fieldtrip or load existing ICA
% By default: baseline correction, scaling, reducing dimensions, whitening, ICA calculation
% Nb of components computation: n_comp = rank(squeeze(data_hilbert_all_nnz.trial{1}) * squeeze(data_hilbert_all_nnz.trial{1})');
if strcmp(compute_ICA, 'yes') == 1
    cfg = [];
    cfg.demean = 'yes';
    cfg.method = 'icasso' %'fastica'; %'icasso' 
    cfg.randomseed = 5;  % any integer, so that the ICA algorithm gives the same results each time it's computed
    cfg.numcomponent = numcomponent;
    cfg.fastica.firstEig = 1;
    cfg.fastica.lastEig = numcomponent;
    % cfg.fastica.interactivePCA = on; % To specify interactively the number of components to retain
    cfg.icasso.mode = 'randinit' % 'both' for bootstrapping the data and randomizing initial conditions
    cfg.icasso.Niter = 100
    ft_ICA = ft_componentanalysis(cfg, data_hilbert_all_nnz);
    tICs_all = ft_ICA.trial{1};
    save([Folder_path '/ICAdecomp/ncomps_' num2str(numcomponent) '/tICs_all_ft_fastica.mat'], 'tICs_all', '-v7.3');
elseif strcmp(compute_ICA, 'no') == 1
    load([Folder_path '/ICAdecomp/ncomps_' num2str(numcomponent) '/icasso_nica15_corr25/tICs_all_ft_fastica.mat']);
end


%% Compute correlation between independent components and raw time course, also make some summary plots
disp('Computing correlation data - tIC')
min_nb_subj_corr = 13;  % Define min number of subjects for which we want to have |r|>=0.3 --> 13 ~ 2/3 of the subjects
maxcorrs = zeros(size(tICs_all,1),1);
mincorrs = zeros(size(tICs_all,1),1);
nbvox_corr = zeros(size(tICs_all,1),1);
nbvox_nocorr = zeros(size(tICs_all,1),1);
nbvox_corr_noglobal = zeros(size(tICs_all,1),2);
vect_cond11 = data_hilbert_all.condition == 11;
vect_cond12 = data_hilbert_all.condition == 12;
vect_cond21 = data_hilbert_all.condition == 21;
vect_cond22 = data_hilbert_all.condition == 22;
vect_condH = logical(vect_cond11 + vect_cond12);
vect_condL = logical(vect_cond21 + vect_cond22);
figure('Color', [1 1 1]);
figure('Color', [1 1 1]);
for t=1:size(tICs_all,1)
    
    fprintf(['  tIC ' num2str(t) ':\n'])
    tIC = tICs_all(t,:);

    % Compute Pearson correlation between tIC and time course of each voxel
    % save result in a Vx1 matrix with N = number of voxels
    corr_tIC_hilbert = zeros(size(data_hilbert_all.trial{1},1),1);  % same with original or modified data
    corr_tIC_hilbert_cond = zeros(size(data_hilbert_all.trial{1},1),2);  
    nb_subj_with_corr = zeros(size(data_hilbert_all.trial{1},1),1);
    nb_subj_no_corr = zeros(size(data_hilbert_all.trial{1},1),1);
    nb_subj_with_corr_noglobal = zeros(size(data_hilbert_all.trial{1},1),2);
    for ngrid = 1:size(data_hilbert_all.trial{1},1)

        % Compute correlation for each voxel
        correlation_coefficients = corrcoef(tIC, data_hilbert_all.trial{1}(ngrid,:));
        corr_tIC_hilbert(ngrid) = correlation_coefficients(2);
        
        % Check how many subjects have the correlation for the current voxel
        if abs(corr_tIC_hilbert(ngrid)) >= threshold_corr
            for isubj = 1:max(vect_subj)
                ivect_subj = vect_subj == isubj;
                correlation_coefficients_subj = corrcoef(tIC(ivect_subj), data_hilbert_all.trial{1}(ngrid,ivect_subj));
                if abs(correlation_coefficients_subj(2)) >= threshold_corr
                    nb_subj_with_corr(ngrid) = nb_subj_with_corr(ngrid)+1;
                else
                    nb_subj_no_corr(ngrid) = nb_subj_no_corr(ngrid)+1;
                end
            end
        else  % Check how many subjects have the correlation for the current voxel when the global correlation is inferior to the reference value
            nb_subj_with_corr(ngrid) = NaN;
            nb_subj_no_corr(ngrid)   = NaN;
            for isubj = 1:max(vect_subj)
                ivect_subj = vect_subj == isubj;
                correlation_coefficients_subj = corrcoef(tIC(ivect_subj), data_hilbert_all.trial{1}(ngrid,ivect_subj));
                if correlation_coefficients_subj(2) >= threshold_corr
                    nb_subj_with_corr_noglobal(ngrid, 1) = nb_subj_with_corr_noglobal(ngrid, 1)+1;
                elseif correlation_coefficients_subj(2) <= - threshold_corr
                    nb_subj_with_corr_noglobal(ngrid, 2) = nb_subj_with_corr_noglobal(ngrid, 2)+1;
                end
            end
        end
        
        % Compute the correlation for each condition separately
        correlation_coefficients = corrcoef(tIC(vect_condH), data_hilbert_all.trial{1}(ngrid,vect_condH));
        corr_tIC_hilbert_cond(ngrid, 1) = correlation_coefficients(2);
        correlation_coefficients = corrcoef(tIC(vect_condL), data_hilbert_all.trial{1}(ngrid,vect_condL));
        corr_tIC_hilbert_cond(ngrid, 2) = correlation_coefficients(2);

    end
    
    % Make histogram of distribution of correlation coefficients for all voxels
    figure(1);
    subplot(5, 5, t);
    hist(corr_tIC_hilbert);
    xlim([-0.65 0.75]);
    ylim([0 2500]);
    title(['tIC ' num2str(t)]);
    
    % Make histogram of distribution of correlation coefficients for all voxels, for each condition
    figure(2);
    subplot(5, 5, t);
    histogram(corr_tIC_hilbert_cond(:, 1), 'facecolor', 'g');
    hold on;
    histogram(corr_tIC_hilbert_cond(:, 2), 'facecolor', 'r');
    line([nanmean(corr_tIC_hilbert_cond(:, 1)) nanmean(corr_tIC_hilbert_cond(:, 1))], [0 850], 'color', 'g', 'linewidth', 2);
    line([nanmean(corr_tIC_hilbert_cond(:, 2)) nanmean(corr_tIC_hilbert_cond(:, 2))], [0 850], 'color', 'r', 'linewidth', 2);
    xlim([-0.65 0.75]);
    ylim([0 850]);
    title(['tIC ' num2str(t)]);
    
    % Select the maximal correlation coefficient obtained across all voxels
    maxcorrs(t)=max(corr_tIC_hilbert);
    mincorrs(t)=min(corr_tIC_hilbert);
    fprintf(['     Maximum correlation: ' num2str(max(corr_tIC_hilbert)) '\n'])
    fprintf(['     Minimum correlation: ' num2str(min(corr_tIC_hilbert)) '\n'])
    
    % Indicate how many voxels have |r|>=0.3 for at least 'min_nb_subj_corr' subjects 
    nbvox_corr(t) = length(find(nb_subj_with_corr>=min_nb_subj_corr));
    nbvox_nocorr(t) = length(find(nb_subj_with_corr<min_nb_subj_corr));
    nbvox_corr_noglobal(t, 1) = length(find(nb_subj_with_corr_noglobal(:, 1)>=min_nb_subj_corr));    
    nbvox_corr_noglobal(t, 2) = length(find(nb_subj_with_corr_noglobal(:, 2)>=min_nb_subj_corr));    
    fprintf(['     Nb of voxels with |r|>=' num2str(threshold_corr) ' for at least ' num2str(min_nb_subj_corr) ' subjects:  ' num2str(nbvox_corr(t)) '\n'])
    fprintf(['     Nb of voxels with |r|>=' num2str(threshold_corr) ' for less than ' num2str(min_nb_subj_corr) ' subjects: ' num2str(nbvox_nocorr(t)) '\n'])
    fprintf(['     Nb of voxels with |r global|<' num2str(threshold_corr) ' but with r>=' num2str(threshold_corr) ' for at least ' num2str(min_nb_subj_corr) ' subjects: ' num2str(nbvox_corr_noglobal(t, 1)) '\n'])
    fprintf(['     Nb of voxels with |r global|<' num2str(threshold_corr) ' but with r<=-' num2str(threshold_corr) ' for at least ' num2str(min_nb_subj_corr) ' subjects: ' num2str(nbvox_corr_noglobal(t, 2)) '\n\n'])
   
end

pct_corr = nbvox_corr * 100 ./ (nbvox_corr+nbvox_nocorr);
figure('color', [1 1 1]);
bar(1:size(tICs_all,1), [nbvox_nocorr,nbvox_corr], 'stacked');
hold on;
text(0.5:size(tICs_all,1)-0.5, nbvox_nocorr+nbvox_corr+30, num2cell(round(pct_corr)));
xlim([0.1 size(tICs_all,1)+1]);
box off;
xlabel('tIC number')
ylabel(['Nb of voxels with |r|>=' num2str(threshold_corr)])
legend(['true for <  ' num2str(min_nb_subj_corr) ' subjects'], ['true for >=' num2str(min_nb_subj_corr) ' subjects']);
title('Comp weight: corr subj by subj');
figure('color', [1 1 1]);
bar(1:size(tICs_all,1), nbvox_corr_noglobal(:, 1));
hold on;
bar(1:size(tICs_all,1), -nbvox_corr_noglobal(:, 2));
xlim([0.1 size(tICs_all,1)+1]);
xlabel('tIC number')
ylabel(['Nb of voxels with |r|<' num2str(threshold_corr) ' but >' num2str(threshold_corr) ' for more than ' num2str(min_nb_subj_corr) ' subjs'])
title(['Strong voxels not leading to general |r|>=' num2str(threshold_corr)])

% Save figures
figure(1); print('-dpng', [Folder_path '/ICAdecomp/ncomps_' num2str(numcomponent) '/Fig_HistComps_Corr_sigma_nnz']);
figure(2); print('-dpng', [Folder_path '/ICAdecomp/ncomps_' num2str(numcomponent) '/Fig_HistComps_Corr_Conds_sigma_nnz']);
figure(3); print('-dpng', [Folder_path '/ICAdecomp/ncomps_' num2str(numcomponent) '/Fig_Weight_Comps_sigma_nnz']);
figure(4); print('-dpng', [Folder_path '/ICAdecomp/ncomps_' num2str(numcomponent) '/Fig_Weight2_Comps_sigma_nnz']);


%% Plot timecourse of average of voxels
nnz_vox = find(data_hilbert_all.trial{1}(:, 1) ~= 0);
avg_vox = mean(data_hilbert_all.trial{1}(nnz_vox, :));
fig_avg_vox = figure('color', [1 1 1]);
plot(avg_vox);
xlim([0 length(avg_vox)]);
figure(fig_avg_vox); print('-dpng', [Folder_path '/ICAdecomp/ncomps_' num2str(numcomponent) '/Fig_TimeCourse_avgVox']);


%% Plot sources of CORRELATIONS of component of interest
template_grid = load([path2ft '/template/sourcemodel/standard_sourcemodel3d10mm.mat']);

% Load template T1 and reslice
templatedir  = [path2ft '/external/spm8/templates'];
template = ft_read_mri(fullfile(templatedir, 'T1.nii'));
template_resliced = ft_volumereslice([], template);

% Select a tIC to compute correlation values across all voxels
ntIC = 15;
threshold_corr = 0.25

tIC = tICs_all(ntIC,:);
corr_tIC_hilbert = zeros(size(data_hilbert_all.trial{1},1),1);
nb_subj_with_corr = zeros(size(data_hilbert_all.trial{1},1),1);
for ngrid = 1:size(data_hilbert_all.trial{1},1)
    
    if data_hilbert_all.trial{1}(ngrid, 1) == 0
        corr_tIC_hilbert(ngrid) = 0;
    else
        correlation_coefficiencts = corrcoef(tIC, data_hilbert_all.trial{1}(ngrid,:));
        corr_tIC_hilbert(ngrid) = correlation_coefficiencts(2);
    end
    
    % Test on nb of subjects having this correlation
    if abs(corr_tIC_hilbert(ngrid)) >= threshold_corr
        for isubj = 1:max(vect_subj)
            ivect_subj = vect_subj == isubj;
            correlation_coefficients_subj = corrcoef(tIC(ivect_subj), data_hilbert_all.trial{1}(ngrid,ivect_subj));
            if abs(correlation_coefficients_subj(2)) >= threshold_corr
                nb_subj_with_corr(ngrid) = nb_subj_with_corr(ngrid)+1;
            end
        end
    else
        nb_subj_with_corr(ngrid) = NaN;
    end
    
end

% Mask correlation values too low and not consistent among subjects (at least 13 subjects showing |r|>=0.3).
%%% Warning: r<-0.3 && r>0.3
corr_tIC_hilbert(abs(corr_tIC_hilbert) < threshold_corr) = 0;  % mask voxels with low correlation values
corr_tIC_hilbert(isnan(corr_tIC_hilbert)) = 0;  % NaN values are when grid point is outside the brain, set them to zero too

% Prepare warning message on command window, to indicate if positive or negative correlation (displayed after figure plotting to be visible on command window)
neg_pts = length(corr_tIC_hilbert(corr_tIC_hilbert<0));
pos_pts = length(corr_tIC_hilbert(corr_tIC_hilbert>0));
if neg_pts > 0 && pos_pts == 0
    corr_tIC_hilbert = corr_tIC_hilbert*-1;
end

% Distinguish voxels not consistent between subjects
corr_tIC_hilbert(nb_subj_with_corr<min_nb_subj_corr) = corr_tIC_hilbert(nb_subj_with_corr<min_nb_subj_corr)*-1;  
    
% Create a fieldtrip-like data structure with the correlation values
correlation_plot = [];
correlation_plot.dim = template_grid.sourcemodel.dim;
correlation_plot.inside = template_grid.sourcemodel.inside;
correlation_plot.pos = template_grid.sourcemodel.pos;
correlation_plot.method = 'average';
correlation_plot.time = data_hilbert_all.time;
correlation_plot.avg.pow = corr_tIC_hilbert; % here are the correltaion values for each grid point

% Interpolation: Project computed sources on standard MRI brain
cfg            = [];
cfg.parameter = 'avg.pow';
correlation_plot_intp  = ft_sourceinterpolate(cfg, correlation_plot, template_resliced);

% Visualize the result on 3 MRI views
%%% Consistent voxels between subjects: warm colors
%%% Inconsistent voxels between subjects: cold colors
cfg              = [];
cfg.method       = 'ortho';
cfg.funparameter = 'avg.pow';
cfg.location = 'min';
cfg.funcolorlim  = 'maxabs';
ft_sourceplot(cfg, correlation_plot_intp);
warning(['This components has: ' num2str(pos_pts) ' positive correlation voxels and ' num2str(neg_pts) ' negative correlation voxels.']);
            

%% Plot particular timecourses

% %%% Plot timecourse of a component of interest
% ntIC = 11;  % define component of interest
% figure('color', [1 1 1]);
% plot(tICs_all(ntIC, :));
% 
% %%% Plot timecourse of the data at a voxel of interest
% nvox = 830;  % define voxel of interest
% figure('color', [1 1 1]);
% plot(data_hilbert_all.trial{1}(nvox, :))
% 
% %%% Plot with a color for each subject
% close all;
% ntIC = 20;
% nvox = [3260 85 130];
% for i=1:length(nvox)+1
%     figure('color', [1 1 1]);
% end
% for isubj = 1:max(vect_subj)
%     ivect_subj = vect_subj == isubj;
%     ibegin = min(find(ivect_subj==1));
%     iend   = max(find(ivect_subj==1));
%     col_subj = [rand rand rand];
%     figure(1);
%     plot(ibegin:iend, tICs_all(ntIC, ivect_subj), 'color', col_subj);
%     hold on;
%     for ivox = 1:length(nvox)
%         figure(ivox+1);
%         plot(ibegin:iend, data_hilbert_all.trial{1}(nvox(ivox), ivect_subj), 'color', col_subj);
%         hold on;
%         title(['Vox' num2str(nvox(ivox))]);
%         xlim([0 iend]);
%     end
% end


%% Plot for one component, the distribution of r across voxels for each subject
% ntIC = 20;
% tIC = tICs_all(ntIC,:);
% figure('color', [1 1 1]);
% for isubj = 1:max(vect_subj)
%     r_subj = zeros(size(data_hilbert_all.trial{1},1), 1);
%     for ngrid = 1:size(data_hilbert_all.trial{1},1)
%         ivect_subj = vect_subj == isubj;
%         correlation_coefficients_subj = corrcoef(tIC(ivect_subj), data_hilbert_all.trial{1}(ngrid,ivect_subj));
%         r_subj(ngrid) = correlation_coefficients_subj(2);
%     end
%     subplot(4, 5, isubj);
%     hist(r_subj);
%     %xlim([-1 1]);
%     ylim([0 4000]);
%     title(['S ' num2str(isubj)]);
% end

   


