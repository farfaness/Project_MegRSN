function [ source_timecourses1cm ] = batch_forwardmodel1cm_beamforming_single_subject(subjects)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This script:
%%%     1) concatenate the alpha band epoched data for the two conditions dataset (2*2runs)
%%%     2) compute a head and grid model
%%%     3) run a LCMV beamforming
%%%
%%%  donne for a particular frequency band, and for each subject on the magnometers channels
%%%  start by loading a preprocessed MRI and a head model for each subj (computed by Batch_MRI_preprocessing_and_head_model)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic

%%%% Indicate the frequency band
name_band_freq = 'theta'
isubj = subjects ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Individual MRI Preprocessing (units conversion and reslicing)

eval(['subject' int2str(isubj)])
Outputdir = eval('subjectdata.subjectdir')
    
     template_grid = load('/lena13/home_users/users/hazem/fieldtrip-20161231/template/sourcemodel/standard_sourcemodel3d10mm.mat');
 
     % Load the sensors structure from the .grad of the first session and first run
     if i == 1 || i == 6 || i == 10 || i == 11 || i == 13 || i == 15 || i == 16 || i == 20 || i == 22;
            sensors_distances_filename = [subjectdata.subjectdir filesep eval('subjectdata.datadir1') filesep 'run01_tsss_realign_ica_eye_corrected_grad_ica_cardio_corrected_grad_epoched_alpha.mat']
            else sensors_distances_filename = [subjectdata.subjectdir filesep eval('subjectdata.datadir1') filesep 'run01_tsss_ica_eye_corrected_grad_ica_cardio_corrected_grad_epoched_alpha.mat']
     end
     load(sensors_distances_filename, 'data_epoched_filtered')
     sensors_distances = data_epoched_filtered.grad 
     clear data_epoched_filtered 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
%%%% Load the realigned and resliced meters MRI, and the head model

load([ eval('Outputdir') filesep 'Anatomy/mri_realigned_resliced_meters.mat'], 'mri_realigned_resliced_meters')
load ([eval('Outputdir') filesep 'Anatomy/headmodel.mat'],'headmodel')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Create the forward model
      
       %Individual-sourcemodel based on a template (warpening)
       %i.e. Make the subject specific grid, using a common template grid
       %(see http://www.fieldtriptoolbox.org/tutorial/sourcemodel#subject-specific_grids_that_are_equivalent_across_subjects_in_normalized_space)
       cfg                 = [];
       cfg.grid.warpmni    = 'yes';
       cfg.grid.template   = template_grid.sourcemodel; %(template en mm, mni convention: http://www.fieldtriptoolbox.org/faq/how_are_the_different_head_and_mri_coordinate_systems_defined) 
       cfg.grid.nonlinear  = 'yes'; % use non-linear normalization
       cfg.mri             = mri_realigned_resliced_meters    % MRI structure, containing the individual anatomy (en meters, mais pas de prblm convertion to spm syst auto)
       cfg.grid.unit       ='mm'; %(by default it takes cm)
 %     cfg.grid.resolution =  % use a 3-D grid with ... resolution, If both are defined cfg.grid.template prevails
       sourcemodel1cm         = ft_prepare_sourcemodel(cfg);
       sourcemodel1cm         = ft_convert_units(sourcemodel1cm, 'm');
 
       save([eval('Outputdir') filesep 'Anatomy/sourcemodel1cm.mat'],'sourcemodel1cm')
       
         % !! make a figure of the single subject headmodel, and grid positions
      figure; hold on;
      ft_plot_vol(headmodel, 'edgecolor', 'none', 'facealpha', 0.4);
      ft_plot_mesh(sourcemodel.pos(sourcemodel.inside,:));
       
       sensors_distances_meters = ft_convert_units(sensors_distances, 'm');
       
       % Create the grid model (Leadfield)
       cfg                 = [];
       cfg.grad            = sensors_distances_meters % structure with sensor definition, see FT_DATATYPE_SENS
       cfg.grid            = sourcemodel1cm;  %(in meters)
       cfg.headmodel       = headmodel;  %(in meters)
       cfg.normalize       = 'yes';
       cfg.reducerank      = 2;
       cfg.grid.unit       = 'm' %'m';
       cfg.channel         = {'megmag'} ; % or {'meggrad'} %Just the magnetometers (Then everything has to be convert to meters, see: https://mailman.science.ru.nl/pipermail/fieldtrip/2017-May/011539.html   or 'meggrad' /  'megplanar' )
       [grid1cm] = ft_prepare_leadfield(cfg);
 
       save([Outputdir filesep 'Anatomy/grid1cm.mat'],'grid1cm')
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%% Concatenate the epoched data for each subject and compute LCMV beamformer

    % CONCATENATE THE EPOCHED DATA
    
    for S = [1 2]  %Pour les 2 sessions
        datadir =['subjectdata.datadir' int2str(S)]
        for R = [1 2]  %Pour les 2 runs
            if isubj == 1 || isubj == 6 || isubj == 10 || isubj == 11 || isubj == 13 || isubj == 15 || isubj == 16 || isubj == 20 || isubj == 22;
                input_filename = [subjectdata.subjectdir filesep eval(datadir) filesep 'run0' int2str(R) '_tsss_realign_ica_eye_corrected_grad_ica_cardio_corrected_grad' '_epoched_' name_band_freq '.mat']
            else input_filename = [subjectdata.subjectdir filesep eval(datadir) filesep 'run0' int2str(R) '_tsss_ica_eye_corrected_grad_ica_cardio_corrected_grad' '_epoched_' name_band_freq '.mat']
            end
            condition = eval(['subjectdata.condition' int2str(S) int2str(R)])
            load(input_filename)
            
            if strcmp(condition,'High_contact') & R == 1
                data_epoched_HC_R1 = data_epoched_filtered
                
            elseif strcmp(condition,'High_contact') & R == 2
                data_epoched_HC_R2 = data_epoched_filtered
                
            elseif strcmp(condition,'Low_contact') & R == 1
                data_epoched_LC_R1 = data_epoched_filtered
                
            elseif strcmp(condition,'Low_contact') & R == 2
                data_epoched_LC_R2 = data_epoched_filtered
                
            end
            clear data_epoched_filtered
        end
    end
    
    %Concatenate the 2 conditions (2*2runs)
    cfg = [];
    data_epoched_concat = ft_appenddata (cfg, data_epoched_HC_R1, data_epoched_HC_R2, data_epoched_LC_R1, data_epoched_LC_R2);
    
    % Add conditions informations to the concatenated data
    n1 = length(data_epoched_HC_R1.trial(1, :))
    n2 = length(data_epoched_HC_R2.trial(1, :))
    n3 = length(data_epoched_LC_R1.trial(1, :))
    n4 = length(data_epoched_LC_R2.trial(1, :))
    
    data_epoched_concat.condition = [1: (n1 + n2 + n3+ n4)]
    data_epoched_concat.condition (1, (1:n1)) = 11
    data_epoched_concat.condition (1, ((1+n1):(n1+n2))) = 12
    data_epoched_concat.condition (1, ((1+n1+n2):(n1+n2+n3))) = 21
    data_epoched_concat.condition (1, ((1+n1+n2+n3):(n1+n2+n3+n4))) = 22
    
    clear data_epoched_HC_R1 data_epoched_HC_R2 data_epoched_LC_R1 data_epoched_LC_R2;
    %----------------------------------------------------------------
    % CONCATENATE THE CONTINOUS DATA
    
    for S = [1 2]  %Pour les 2 sessions
        datadir =['subjectdata.datadir' int2str(S)]
        for R = [1 2]  %Pour les 2 runs
            if isubj == 1 || isubj == 6 || isubj == 10 || isubj == 11 || isubj == 13 || isubj == 15 || isubj == 16 || isubj == 20 || isubj == 22;
                input_filename = [subjectdata.subjectdir filesep eval(datadir) filesep 'run0' int2str(R) '_tsss_realign_ica_eye_corrected_grad_ica_cardio_corrected_grad' '_continuous_' name_band_freq '.mat']
            else input_filename = [subjectdata.subjectdir filesep eval(datadir) filesep 'run0' int2str(R) '_tsss_ica_eye_corrected_grad_ica_cardio_corrected_grad' '_continuous_' name_band_freq '.mat']
            end
            condition = eval(['subjectdata.condition' int2str(S) int2str(R)])
            load(input_filename)
            
            if strcmp(condition,'High_contact') & R == 1
                data_continuous_HC_R1 = data_continuous_filtered
                
            elseif strcmp(condition,'High_contact') & R == 2
                data_continuous_HC_R2 = data_continuous_filtered
                
            elseif strcmp(condition,'Low_contact') & R == 1
                data_continuous_LC_R1 = data_continuous_filtered
                
            elseif strcmp(condition,'Low_contact') & R == 2
                data_continuous_LC_R2 = data_continuous_filtered
                
            end
            clear data_continuous_filtered
        end
    end
    
    %Concatenate the 2 conditions (2*2runs)
    cfg = [];
    %       data_continuous_concat       = ft_appenddata (cfg, data_continuous_HC_R1, data_continuous_HC_R2, data_continuous_LC_R1, data_continuous_LC_R2);
    data_continuous_concat       = data_continuous_HC_R1;
    data_continuous_concat.trial = {[data_continuous_HC_R1.trial{1} data_continuous_HC_R2.trial{1} data_continuous_LC_R1.trial{1} data_continuous_LC_R2.trial{1}]};
    data_continuous_concat.time  = {[0:1:1*((size(data_continuous_HC_R1.time{1},2) + size(data_continuous_HC_R2.time{1},2) + size(data_continuous_LC_R1.time{1},2) + size(data_continuous_LC_R2.time{1},2)) -1)]};
    
    % Add conditions informations to the concatenated data
    n1 = length(data_continuous_HC_R1.time{1, 1})
    n2 = length(data_continuous_HC_R2.time{1, 1})
    n3 = length(data_continuous_LC_R1.time{1, 1})
    n4 = length(data_continuous_LC_R2.time{1, 1})
    
    data_continuous_concat.condition = [1: (n1 + n2 + n3+ n4)]
    data_continuous_concat.condition (1, (1:n1)) = 11
    data_continuous_concat.condition (1, ((1+n1):(n1+n2))) = 12
    data_continuous_concat.condition (1, ((1+n1+n2):(n1+n2+n3))) = 21
    data_continuous_concat.condition (1, ((1+n1+n2+n3):(n1+n2+n3+n4))) = 22
    
    clear data_continuous_HC_R1 data_continuous_HC_R2 data_continuous_LC_R1 data_continuous_LC_R2;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %% COMPUTE THE BEAMFORMER
    
%     %Load the grid and the headmodel
%     load([eval('Outputdir') filesep 'Anatomy/headmodel.mat'])
%     grid1cm = load([eval('Outputdir') filesep 'Anatomy/grid.mat'])
    
    
%LCMV beamformer : calculate sources for the concatenated datasets (Apply a lcmv beamformer and extract the source timecourses)
% compute covariance matrix : average all the trials in a data structure and also estimates the covariance
% !!!! For a correct covariance estimation it is important that you used the cfg.demean = 'yes' option when the function ft_preprocessing was applied
fprintf('\n###############\nComputing covariance matrix...\n\n')
cfg                  = [];
cfg.channel          = {'megmag'} 
cfg.covariance       = 'yes';
cfg.covariancewindow = 'all';  %'prestim', 'poststim', 'all' or [begin end] (default = 'all')
cfg.vartrllength     = 2; %accept variable length trials, use all available trials the available samples in every trial will be used for the average and covariance computation. Missing values are replaced by NaN and are not included in the computation.
timelock       = ft_timelockanalysis(cfg, data_epoched_concat);

% % extract MEG channels
chansel = ft_channelselection('megmag', data_continuous_concat.label); % find MEG sensor names
chansel = match_str(data_continuous_concat.label, chansel);         % find MEG sensor indices

% initialize source time course structure
source_timecourses1cm   = [];
source_timecourses1cm.pos     = grid1cm.pos;
source_timecourses1cm.trial   = [];
source_timecourses1cm.mom     = []; % here the strength of the sources in the three directions will be stored
source_timecourses1cm.filter  = []; % here the beamformer filter for the respective grid will be stored
source_timecourses1cm.label   = {'brain'} %{'cortex'};

% downsample continuous data to 100 Hz for faster computation
fprintf('\n###############\nDownsampling data to 100 Hz...\n\n')
cfg         = [];
cfg.resamplefs = 100;
cfg.detrend    = 'no';
cfg.demean     = 'no';
data_continuous_concat_down = ft_resampledata(cfg, data_continuous_concat);

source_timecourses1cm.time = data_continuous_concat_down.time;

% iterate over grids in the volume, compute filter, mom and trial
fprintf('\n###############\nExtracting time courses for each grid...\n\n')
for ngrid=1:size(grid1cm.pos, 1) 
    
    % check if voxel inside the brain
    if grid1cm.inside(ngrid) == 1
        
        fprintf(['\n###############\nScanning grid ' num2str(ngrid) ' of ' num2str(size(grid1cm.pos, 1)) ' ('  num2str(length(find(grid1cm.inside ==1))) ' inside the volume)...\n\n'])
        % calculate the beamformer filter
        cfg                   = [];
        cfg.method            = 'lcmv';
        cfg.grad              = data_continuous_concat.grad ;
        cfg.headmodel         = headmodel;
        cfg.grid.pos          = grid1cm.pos(ngrid, :);
        cfg.grid.inside       = 1:size(cfg.grid.pos, 1);
        cfg.grid.outside      = [];
        cfg.lcmv.keepfilter   = 'yes';
        source_idx            = ft_sourceanalysis(cfg, timelock);
        beamformer            = source_idx.avg.filter{1};
        
        source_timecourses1cm.filter = [source_timecourses1cm.filter; {beamformer}];
        
        % compute mom (source power in each direction)
        pow_data = [];
        pow_data.label = {'pow_x', 'pow_y', 'pow_z'};
        pow_data.time = data_continuous_concat_down.time;
        for isubj=1:length(data_continuous_concat_down.trial) 
            pow_data.trial{isubj} = beamformer * data_continuous_concat_down.trial{isubj}(chansel,:);
        end
        
        source_timecourses1cm.mom = [source_timecourses1cm.mom; {pow_data.trial{isubj}}];

        % project the source to its strongest orientation, i.e. the 
        % direction that explains most of the source variance. This is done
        % by single-value decomposition that is equivalent to PCA
        mom_voxel = cat(2, pow_data.trial{:});
        [u, s, v] = svd(mom_voxel, 'econ');

        % this is equal to the first column of matrix V, apart from the scaling with s(1,1)
        timeseries_voxel = u(:,1)' * mom_voxel;
        source_timecourses1cm.trial = [source_timecourses1cm.trial; timeseries_voxel];
    
    else
        source_timecourses1cm.trial = [source_timecourses1cm.trial; zeros(1,length(data_continuous_concat_down.time{1}))];
    end
    
end

source_timecourses1cm.trial = {source_timecourses1cm.trial};
save([subjectdata.subjectdir filesep 'HC_and_LC_' name_band_freq '_source_timecourses1cm.mat'], 'source_timecourses1cm', '-v7.3')

clearvars -except source_timecourses1cm

e = toc;
e = num2str(e);
fprintf(['###################   The batch Forward model 1 cm and Beamforming took ' e ' seconds #######################']);

end


