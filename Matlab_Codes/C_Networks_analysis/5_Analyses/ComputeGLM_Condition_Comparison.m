%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute GLM
%%%
%%% different methods for analysis possibles:
%%% 'completeGLM': timecourse = bcomp1*comp1 + bcomp2*comp2 + ... + bcompN*compN
%%% 'linreg': timecourse = bcompi*compi
%%% 'luckhoo': compi = cst + betacond1*indexcond1 + e
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc

%% Init
freq_band = 'theta_b'; % 'alpha' or 'beta'
method_analysis = 'linreg';
ncompsICA = 15;  % 25 or 18

%%% Methods for analysis, choose among:
% 'completeGLM': timecourse = bcomp1*comp1 + bcomp2*comp2 + ... + bcompN*compN
% 'linreg': timecourse = bcompi*compi
% 'luckhoo': compi = cst + betacond1*indexcond1 + e


%% Set path
Folder_path = ['/pclxserver/raid11/data/SELFI/Networks_Group_comparison/' freq_band];
path2ft = '/lena13/home_users/users/hazem/fieldtrip-20161231'
addpath(path2ft);
ft_defaults


%% ICA components
% Load data (computed in step 4_ICA)
if ncompsICA == 18
    disp('###### Loading ICA components: n=18, fastica method');
    load([Folder_path '/ICAdecomp/ncomps_18/tICs_all_ft_fastica.mat']);
    Folder_path_results = [Folder_path '/ICAdecomp/ncomps_18'];
elseif ncompsICA == 25
    disp('###### Loading ICA components: n=25, fastica method');
    load([Folder_path '/ICAdecomp/tICs_all_ft_fastica.mat']);
    Folder_path_results = [Folder_path '/ICAdecomp'];
elseif ncompsICA == 15
    disp('###### Loading ICA components: n=15, icasso method');
    load([Folder_path '/ICAdecomp/ncomps_15/icasso_nica15_corr25/tICs_all_ft_fastica.mat']);
    Folder_path_results = [Folder_path '/ICAdecomp/ncomps_15/icasso_nica15_corr25'];
end


% load('/pclxserver/raid11/data/SELFI/Networks_Group_comparison/alpha/ICAdecomp/ncomps_15/icasso_nica15_corr25/tICs_all_ft_fastica.mat')


% Select components of interest to include in the GLM
%%% Warning: procedure to be improved
if strcmp(freq_band, 'alpha') == 1  
    if ncompsICA == 18
        comps_interest = [1 3 4 6 7 8 10 11 12 13 15 18];
    elseif ncompsICA == 25
        comps_interest = [4 17 24 25];%[2, 8, 9, 11, 12, 16, 17, 24, 25];  % >15% of voxels had at least 13 subj with |rsubj|>0.3
    elseif ncompsICA == 15
        comps_interest = [2:7, 9:11, 13, 15];%1:15;
    end
elseif strcmp(freq_band, 'beta') == 1  
    comps_interest = [12, 18, 22, 23];  % >10% of voxels had at least 13 subj with |rsubj|>0.3
elseif strcmp(freq_band, 'theta_b') == 1  
    comps_interest = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15];  
end


%% Load concatenated source data to extract subj and cond information
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
nsubj = max(vect_subj);

% Create vector indicating for each data point to which condition it corresponds to
vect_cond = data_hilbert_all.condition;

% Create consistent label field
for i=1:size(data_hilbert_all.trial{1}, 1)
    data_hilbert_all.label{i} = ['v' num2str(i)];
end


%% ANALYZE with methods 'completeGLM' and 'linreg'
if strcmp(method_analysis, 'completeGLM') == 1 || strcmp(method_analysis, 'linreg') == 1
    betas = [];
    betas.cond1 = zeros(length(comps_interest)+1, size(data_hilbert_all.trial{1}, 1), nsubj);  % +1 because of column for intercept
    betas.cond2 = zeros(length(comps_interest)+1, size(data_hilbert_all.trial{1}, 1), nsubj);  % +1 because of column for intercept
    for isubj = 1:nsubj
        
        for icond = 1:2
            
            % Feedback on screen
            disp(['Running isubj ' num2str(isubj) ', cond ' num2str(icond)]);
            
            % Set coding of conditions
            if icond == 1
                % High contact, blocks 1 and 2
                val_cond = [11 12];
            elseif icond == 2
                % Low contact, blocks 1 and 2
                val_cond = [21 22];
            end
            
            % Find indices of data corresponding to current subject and current condition
            icond_data = find(vect_cond == val_cond(1) | vect_cond == val_cond(2));
            isubj_data = find(vect_subj == isubj);
            idata = intersect(icond_data, isubj_data);
            
            % Build predictor matrix (time x components of interest) and zscore across time points
            predict = zeros(length(idata), length(comps_interest));
            for ipred = 1:length(comps_interest)
                predict(:, ipred) = tICs_all(comps_interest(ipred), idata);
            end
            zpredict = predict;%zscore(predict); % no need of zscoring
            
            % Compute for each voxel inside the brain
            for ivox = 1:size(data_hilbert_all.trial{1}, 1)
                
                % Voxel outside the brain --> NaN
                if data_hilbert_all.trial{1}(ivox, 1) == 0
                    betas.(['cond' num2str(icond)])(:, ivox, isubj) = NaN;
                    
                % Voxel inside the brain --> compute GLM
                else
                    % Build Y vector for the current voxel (source timecourse, for current subj and cond) and zscore across time points
                    Y = data_hilbert_all.trial{1}(ivox, idata);
                    zY = Y;%zscore(Y); % no need of zscoring
                    
                    % METHOD complete GLM with all components
                    if strcmp(method_analysis, 'completeGLM') == 1
                        % Compute GLM and store output beta values (first beta corresponds to constant term)
                        B = glmfit(zpredict, zY, 'normal');     % Same as with Catherine's function: [B, cst, ~] = My_glm(zY', zpredict, 0); B = [cst; B];
                        betas.(['cond' num2str(icond)])(:, ivox, isubj) = B;
                        
                    % METHOD one regression per component
                    elseif strcmp(method_analysis, 'linreg') == 1
                        betas.(['cond' num2str(icond)])(1, ivox, isubj) = NaN;  % this corresponds to the intercept in the complete glm model
                        % Compute linear regression for each component
                        for icomp = 1:size(zpredict, 2)
                            B = regress(zY', zpredict(:, icomp));
                            betas.(['cond' num2str(icond)])(icomp+1, ivox, isubj) = B;
                        end
                        
                    end
                end
                
            end
            
        end
        
    end
    
    %%% Create neighbours for clustering
    disp('Creating neighbouring matrix...');
    cfg = [];
    cfg.method = 'distance';
    cfg.neighbourdist = 0.01;
    cfg.elec.elecpos = data_hilbert_all.pos;
    cfg.elec.chanpos = data_hilbert_all.pos;
    cfg.elec.label   = data_hilbert_all.label;
    neighbours = ft_prepare_neighbours(cfg);
    

    
%%% Compute clustering for each component, with correction for multiple comparisons over components tested
nrands = 5000;
posdistrib_all = zeros(length(comps_interest), nrands);
negdistrib_all = zeros(length(comps_interest), nrands);

    for icomp = 1:length(comps_interest)
        
        % Build fieldtrip-like data structures (grand average on sensors) for clustering
        % icomp+1: because the first component is the intercept in the completeGLM model
        data_cond1 = [];
        data_cond1.individual = permute(squeeze(betas.cond1(icomp+1, :, :)), [1, 3, 2]);
        data_cond1.label = data_hilbert_all.label;
        data_cond1.dimord = 'chan_time_subj';
        data_cond1.time = 1;
        data_cond2 = data_cond1;
        data_cond2.individual = permute(squeeze(betas.cond2(icomp+1, :, :)), [1, 3, 2]);
        
        % Compute clustering
        cfg = [];
        cfg.method = 'montecarlo';
        cfg.statistic = 'depsamplesT';
        cfg.correctm = 'cluster';
        cfg.clusteralpha = 0.025;
        cfg.clusterstatistic = 'maxsum';
        cfg.minnbchan = 0 %3;
        cfg.tail = 0;
        cfg.clustertail = 0;
        cfg.alpha = 0.05;
        cfg.correcttail = 'prob';
        cfg.numrandomization = nrands;
        cfg.neighbours = neighbours;
        cfg.randomseed = 1;  % choose a number to initialize the random seed, so that all clusterings will have the same seed and so the same randomizations
        design = zeros(2, 2*nsubj);
        for i=1:nsubj
            design(1, i) = i;
        end
        for i = 1:nsubj
            design(1, nsubj+i) = i;
        end
        design(2, 1:nsubj) = 1;
        design(2, nsubj+1:2*nsubj) = 2;
        cfg.design = design;
        cfg.uvar = 1;
        cfg.ivar = 2;
        stat_comp.(['comp' num2str(icomp)]) = ft_timelockstatistics(cfg, data_cond1, data_cond2);
        stat_comp.(['comp' num2str(icomp)]).comp_nb = comps_interest(icomp);
        
        % Store all max sum t distributions
        if isfield(stat_comp.(['comp' num2str(icomp)]), 'posdistribution')
            posdistrib_all(icomp, :) = stat_comp.(['comp' num2str(icomp)]).posdistribution;
        end
        if isfield(stat_comp.(['comp' num2str(icomp)]), 'negdistribution')
            negdistrib_all(icomp, :) = stat_comp.(['comp' num2str(icomp)]).negdistribution;
        end
        
    end
    
    % Build new max sum t distributions
    if size(posdistrib_all, 1) > 1
        posdistrib_new = max(posdistrib_all);
        negdistrib_new = min(negdistrib_all);
    else
        posdistrib_new = posdistrib_all;
        negdistrib_new = negdistrib_all;
    end
    
    % Compute new Montecarlo p-values, corrected for multiple comparisons over all the components tested
    for icomp = 1:length(comps_interest)
        if isfield(stat_comp.(['comp' num2str(icomp)]), 'posdistribution')
            if ~isempty(stat_comp.(['comp' num2str(icomp)]).posclusters)
                for iclust = 1:length(stat_comp.(['comp' num2str(icomp)]).posclusters)
                    istat = stat_comp.(['comp' num2str(icomp)]).posclusters(iclust).clusterstat;
                    n_larger = length(find(posdistrib_new >= istat));
                    stat_comp.(['comp' num2str(icomp)]).posclusters(iclust).prob_corr = (n_larger+1) / (nrands+1) * 2;  % see formula in: clusterstat.m, raw 418 (+1: because includes original t-test, *2: because two-tail)
                end
            end
        end
        if isfield(stat_comp.(['comp' num2str(icomp)]), 'negdistribution')
            if ~isempty(stat_comp.(['comp' num2str(icomp)]).negclusters)
                for iclust = 1:length(stat_comp.(['comp' num2str(icomp)]).negclusters)
                    istat = stat_comp.(['comp' num2str(icomp)]).negclusters(iclust).clusterstat;
                    n_larger = length(find(negdistrib_new <= istat));
                    stat_comp.(['comp' num2str(icomp)]).negclusters(iclust).prob_corr = (n_larger+1) / (nrands+1) * 2;
                end
            end
        end
    end
    
    % Save statistics results
    save([Folder_path_results '/stat_comp_5000rands_cleancorr_clusteralpha025_minchan0.mat'], 'stat_comp', 'posdistrib_new', 'negdistrib_new', '-v7.3');
    
    % Display maps of p-values (from clustering) for each component
    % Option with standard MRI (blurry images)
    template_grid = load([path2ft '/template/sourcemodel/standard_sourcemodel3d10mm.mat']);
    templatedir  = [path2ft '/external/spm8/templates'];
    template = ft_read_mri(fullfile(templatedir, 'T1.nii'));
    template_resliced = ft_volumereslice([], template);
    indiv_mri = 1;
    template_mri = ft_read_mri([Folder_path_results '/Masks_nii/Yeo_JNeurophysiol11_MNI152/FSL_MNI152_FreeSurferConformed_1mm.nii']);
    template_mri_resliced = ft_volumereslice([], template_mri);
    template_mri_grey = ft_read_mri([Folder_path_results '/Masks_nii/Yeo_JNeurophysiol11_MNI152/Template_Segmented_SPM12/c1FSL_MNI152_FreeSurferConformed_1mm.nii']);
    template_mri_grey_resliced = ft_volumereslice([], template_mri_grey);
    
    % Option with MRI from one subject (clearer images and allows to apply mask on image)
%     template_grid = load('../../../selfi01_s01/Anatomy/grid1cm.mat');
%     template_grid.sourcemodel = template_grid.grid1cm;
%     template_grid = rmfield(template_grid, 'grid1cm');
%     template_resliced = load('../../../selfi01_s01/Anatomy/mri_realigned_resliced_meters.mat');
%     template_resliced = template_resliced.mri_realigned_resliced_meters;
%     load('/pclxserver/raid11/data/SELFI/selfi01_s01/Anatomy/segmentedmri.mat');
%     template_resliced.anatomy = template_resliced.anatomy .* segmentedmri.template_resliced;
%     ind_mri = 1;

    for icomp = 1:length(comps_interest)
        if ~isempty(find(stat_comp.(['comp' num2str(icomp)]).prob < 0.15))  % Plot maps having significant results
            
            % Feedback on screen
            disp('############################################');
            disp(['Comp ' num2str(comps_interest(icomp)) ':']);
            disp(['   Min uncorr pval = ' num2str(min(stat_comp.(['comp' num2str(icomp)]).prob))]);
            
            % Get p-values for this component map
            plot_tIC = stat_comp.(['comp' num2str(icomp)]).stat; % With voxel by voxel t-test: plot_tIC = map_pval(icomp, :);
            max_pos = max(stat_comp.(['comp' num2str(icomp)]).stat);
            max_neg = abs(min(stat_comp.(['comp' num2str(icomp)]).stat));
            if max_pos > max_neg
                plot_tIC(stat_comp.(['comp' num2str(icomp)]).posclusterslabelmat ~= 1) = 0; % mask pvalues not in the cluster of interest
                str_clust = 'pos cluster';
                p_val = stat_comp.(['comp' num2str(icomp)]).posclusters(1).prob_corr;
                disp(['   Corrected across comps, first cluster pval = ' num2str(p_val)]);
            else
                plot_tIC(stat_comp.(['comp' num2str(icomp)]).negclusterslabelmat ~= 1) = 0; % mask pvalues not in the cluster of interest
                str_clust = 'neg cluster';
                p_val = stat_comp.(['comp' num2str(icomp)]).negclusters(1).prob_corr;
                disp(['   Corrected across comps, first cluster pval = ' num2str(p_val)]);
            end
            disp('############################################');
            
            % Create a fieldtrip-like data structure with the pvalues
            pval_plot = [];
            pval_plot.dim = template_grid.sourcemodel.dim;
            pval_plot.inside = template_grid.sourcemodel.inside;
            pval_plot.pos = template_grid.sourcemodel.pos;
            pval_plot.method = 'average';
            pval_plot.time{1} = data_hilbert_all.time;
            pval_plot.pow = plot_tIC; % here are the p-values for each grid point
            
            % Interpolation: Project on standard MRI brain
            cfg             = [];
            cfg.parameter   = 'pow';
            pval_plot_intp  = ft_sourceinterpolate(cfg, pval_plot, template_resliced);
            
            % Mask non significant values (t-val > 2)
            if indiv_mri == 1
                pval_plot_intp.mask = pval_plot_intp.pow > 2.1 & ~isnan(pval_plot_intp.pow);
                mask_mat = reshape(pval_plot_intp.mask, [256 256 256]) ;
                mask_mat = mask_mat .* segmentedmri.brain ;
                pval_plot_intp.mask = reshape(mask_mat,[256*256*256 1]) ;
            end
            
            % Visualize the result on 3 MRI views
            cfg              = [];
            cfg.method       = 'ortho';
            cfg.funparameter = 'pow';
            cfg.maskparameter = 'mask';
            ft_sourceplot(cfg, pval_plot_intp);
            title(['Comp ' num2str(comps_interest(icomp)) ', ' str_clust ': p=' num2str(p_val)]);
            

            
        else
            
            disp('    ');
            disp('    ');
            disp('############################################');
            disp(['Comp ' num2str(comps_interest(icomp)) ': no significant clusters.']);
            disp(['   Min uncorr pval = ' num2str(min(stat_comp.(['comp' num2str(icomp)]).prob))]);
            disp('############################################');
            disp('    ');
            disp('    ');
            
        end
    end
end



%% ANALYZE with method 'luckhoo'
if strcmp(method_analysis, 'luckhoo') == 1
    
    % Init var to store results
    betas = zeros(nsubj, length(comps_interest), 2);
    betas_envHilb = zeros(nsubj, 2);
    
    for isubj = 1:nsubj
        
        % Find indices corresponding to data of the current subject
        isubj_data = find(vect_subj == isubj);
        
        % Find indices corresponding to high and low conditions
        cond1 = find(vect_cond(isubj_data) == 11 | vect_cond(isubj_data) == 12);
        cond2 = find(vect_cond(isubj_data) == 21 | vect_cond(isubj_data) == 22);
        
        % Compute regression for each component:
        % --> component timecourse = cst + betacond1*flagcond1 + e
        %     cst: beta for condition 2, betacond1: beta for condition 1
        for icomp = 1:length(comps_interest)
            % Analyze ICAcomponents
            data_comp = tICs_all(comps_interest(icomp), isubj_data);
            zdata_comp = zscore(data_comp);
            flag_cond1 = zeros(length(isubj_data), 1);
            flag_cond1(cond1) = 1;
            B = glmfit(flag_cond1, zdata_comp, 'normal');  % 1st: cst (beta cond2), 2nd: beta (cond1)
            betas(isubj, icomp, :) = B;
        end
        
        % Analyze timecourse
        map_betas = zeros(2, nsubj, size(data_hilbert_all.trial{1}, 1));  
        for ivox = 1:size(data_hilbert_all.trial{1}, 1)
            
            % Only compute for voxels inside the brain
            % Voxels outside the brain --> 0
            if data_hilbert_all.trial{1}(ivox, 1) ~= 0
                data_envHilb = data_hilbert_all.trial{1}(ivox, isubj_data);
                zdata_envHilb = zscore(data_envHilb);
                flag_cond1 = zeros(length(isubj_data), 1);
                flag_cond1(cond1) = 1;
                B = glmfit(flag_cond1, zdata_envHilb, 'normal');  % 1st: cst (beta cond2), 2nd: beta (cond1)
                map_betas(:, isubj, ivox) = B;
            end
            
        end
        
    end
    
    % T-test betas on ICs
    disp('  ');
    disp('Results of analysis (Luckhoo''s)');
    for icomp = 1:length(comps_interest)
        
        % Compute t-test
        [~, p] = ttest(squeeze(betas(:, icomp, 1)), squeeze(betas(:, icomp, 2)));
        
        % Feedback on command window
        disp(['- Comp ' num2str(comps_interest(icomp)) ': p=' num2str(p)]);
        disp(['    Mean beta HighContact: ' num2str(mean(betas(:, icomp, 2)))]);
        disp(['    Mean beta LowContact : ' num2str(mean(betas(:, icomp, 1)))]);
        disp('  ');
        
    end
    
    % T-test betas on timecourses
    map_pvalues = zeros(size(map_betas, 3), 1);
    for ivox = 1:size(data_hilbert_all.trial{1}, 1)
        
        % Only compute for voxels inside the brain
        % Voxels outside the brain --> p=1
        if data_hilbert_all.trial{1}(ivox, 1) ~= 0
            [~, p] = ttest(map_betas(1, :, ivox), map_betas(2, :, ivox));
            map_pvalues(ivox) = p;
        else
            map_pvalues(ivox) = 1;
        end
        
    end
    disp(['Nb of significant voxels:  ' num2str(length(find(map_pvalues <= 0.05)))]);
    
    % Display maps of p-values for timecourse
    template_grid = load([path2ft '/template/sourcemodel/standard_sourcemodel3d10mm.mat']);
    templatedir  = [path2ft '/external/spm8/templates'];
    template = ft_read_mri(fullfile(templatedir, 'T1.nii'));
    template_resliced = ft_volumereslice([], template);
    
    % Create a fieldtrip-like data structure with the pvalues
    pval_plot = [];
    pval_plot.dim = template_grid.sourcemodel.dim;
    pval_plot.inside = template_grid.sourcemodel.inside;
    pval_plot.pos = template_grid.sourcemodel.pos;
    pval_plot.method = 'average';
    pval_plot.time{1} = data_hilbert_all.time;
    pval_plot.avg.pow = map_pvalues; % here are the p-values for each grid point
    
    % Interpolation: Project on standard MRI brain
    cfg             = [];
    cfg.parameter   = 'pow';
    pval_plot_intp  = ft_sourceinterpolate(cfg, pval_plot, template_resliced);
    
    % Mask non significant values
    pval_plot_intp.pow(pval_plot_intp.pow > 0.05) = 0.06;
    
    % Visualize the result on 3 MRI views
    cfg              = [];
    cfg.method       = 'ortho';
    cfg.funparameter = 'pow';
    ft_sourceplot(cfg, pval_plot_intp);
    
end




















