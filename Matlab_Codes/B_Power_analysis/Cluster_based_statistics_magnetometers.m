%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Cluster_based statistics
%%%%
%%%% on  magnetometers channels
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%%%% Compute the contrast Condition (High_contact versus Low_contact) in channels * freq space

clear all

load('/pclxserver/raid11/data/SELFI/scripts/neighbours/neuromag306mag_neighb.mat','neighbours');
load( '/pclxserver/raid11/data/SELFI/Grand_average/Power_ind_HighContact.mat', 'power_high_contact_ind' );
load( '/pclxserver/raid11/data/SELFI/Grand_average/Power_ind_LowContact.mat', 'power_low_contact_ind' );

cfg = []
cfg.channel = 'MEG*1';
%cfg.latency     = [begin end] in seconds or 'all' (default = 'all')
cfg.frequency   = [1 40]; 
%cfg.avgoverchan = 'yes' %or 'no'                   (default = 'no')
%cfg.avgovertime = 'yes' or 'no'                   (default = 'no')
%cfg.avgoverfreq = 'yes' or 'no'                   (default = 'no')
cfg.parameter   = 'powspctrm';

cfg.method      =  'montecarlo';     %different methods for calculating the significance probability and/or critical value
                                     %'montecarlo'    get Monte-Carlo estimates of the significance probabilities and/or critical values from the permutation distribution,
                                     %'analytic'      get significance probabilities and/or critical values from the analytic reference distribution (typically, the sampling distribution under the null hypothesis),
                                     %'stats'         use a parametric test from the MATLAB statistics toolbox,
                                     %'crossvalidate' use crossvalidation to compute predictive performance
                                    
cfg.statistic = 'ft_statfun_depsamplesT'   %'ft_statfun_indepsamplesT'for independent samples T-statistic as a measure to evaluate the effect at the sample level
                                           %'ft_statfun_depsamplesT'for within-UO design, dependent samples T-statistic 
                                           %'ft_statfun_indepsamplesF' to compare more than two experimental conditions
                                           %'ft_statfun_depsamplesFmultivariate' to compare more than two experimental conditions
cfg.correctm = 'cluster';
cfg.clusterthreshold = 'parametric';       %method for single-sample threshold, 'parametric', 'nonparametric_individual', 'nonparametric_common' (default = 'parametric')

cfg.clusteralpha = 0.05;  % alpha level of the sample-specific test statistic that will be used for thresholding
cfg.clustertail = 0;
cfg.clusterstatistic = 'maxsum';           %sum of the sample-specific T-statistics that belong to this cluster. Taking the largest of these cluster-level statistics of the different clusters produces the actual test statistic. 
                                           %'wcm' refers to 'weighted cluster mass', a statistic that combines cluster size and intensity; see Hayasaka & Nichols (2004) NeuroImage
cfg.minnbchan = 2;                         % minimum number of neighborhood channels that is required for a selected sample to be included  in the clustering algorithm (default=0).                                 
cfg.neighbours = neighbours;               %neighbourhood structure, see FT_PREPARE_NEIGHBOURS

cfg.tail = 0;                              % -1, 1 or 0 (default = 0); one-sided or two-sided test
cfg.alpha = 0.025;                         % alpha level of the permutation test
cfg.numrandomization = 1000;              % number of draws from the permutation distribution

%Design matrix: 
subj = 19;

design = zeros(2,2*subj);
for i = 1:subj
  design(1,i) = i;
end
for i = 1:subj
  design(1,subj+i) = i;
end
design(2,1:subj)        = 1;
design(2,subj+1:2*subj) = 2;


cfg.design = design;             % design matrix
cfg.uvar = 1;                    % unit of observation (subjects or trials), this allow paired comparison
cfg.ivar  = 2;                   % number or list with indices indicating the independent variable(s)

stat = ft_freqstatistics(cfg, power_high_contact_ind,  power_low_contact_ind)

%Save_the_output
save('/pclxserver/raid11/data/SELFI/Grand_average/stat_cluster_based_High_vs_Low_contact', 'stat')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot the results

cfg = [];
cfg.alpha  = 0.025;
cfg.parameter = 'stat';
cfg.zlim   = [-4 4];
cfg.layout = '/lena13/home_users/users/hazem/fieldtrip-20161231/template/layout/neuromag306mag.lay';
ft_clusterplot(cfg, stat);


