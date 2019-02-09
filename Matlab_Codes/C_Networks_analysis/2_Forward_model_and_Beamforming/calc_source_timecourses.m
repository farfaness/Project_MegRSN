function [ source_timecourses ] = calc_source_timecourses( data_epoched_concat, data_continuous_concat, grid_tmp, headmodel )
%CALC_SOURCE_TIMECOURSES 
%
%   Calculates the source time courses of the data
%   Computes a LCMV beamformer filter based on the epoched data, and applies
%   this filter to the continous data, yielding the strength of the source
%   in each direction (mom). 
%
%   Then computes the time courses on a voxel-by-voxel basis
%   by using a single-value-decomposition that is equivalent to PCA.
%   
%   Saves the result in a struct containing voxel position and voxel time courses.
%   
%   
%   INPUT
%   data_epoched_concat - the epoched bandpass filtered data
%   data_continuous_concat - the continuous bandpass filtered data
%   grid
%   headmodel
% 
%   OUTPUT
%   source_timecourses - the struct containing source positions (source_timecourses.pos)
%                        and timecourses (source_timecourses.trial)


%% tic

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
source_timecourses   = [];
source_timecourses.pos     = grid_tmp.pos;
source_timecourses.trial   = [];
source_timecourses.mom     = []; % here the strength of the sources in the three directions will be stored
source_timecourses.filter  = []; % here the beamformer filter for the respective grid will be stored
source_timecourses.label   = {'brain'} %{'cortex'};

% downsample continuous data to 100 Hz for faster computation
fprintf('\n###############\nDownsampling data to 100 Hz...\n\n')
cfg         = [];
cfg.resamplefs = 100;
cfg.detrend    = 'no';
cfg.demean     = 'no';
data_continuous_concat_down = ft_resampledata(cfg, data_continuous_concat);

source_timecourses.time = data_continuous_concat_down.time;

% iterate over grids in the volume, compute filter, mom and trial
fprintf('\n###############\nExtracting time courses for each grid...\n\n')
for ngrid=1:size(grid_tmp.pos, 1) 
    
    % check if voxel inside the brain
    if grid_tmp.inside(ngrid) == 1
        
        fprintf(['\n###############\nScanning grid ' num2str(ngrid) ' of ' num2str(size(grid_tmp.pos, 1)) ' ('  num2str(length(find(grid_tmp.inside ==1))) ' inside the volume)...\n\n'])
        % calculate the beamformer filter
        cfg                   = [];
        cfg.method            = 'lcmv';
        cfg.grad          = data_continuous_concat.grad ;
        cfg.headmodel         = headmodel;
        cfg.grid.pos          = grid_tmp.pos(ngrid, :);
        cfg.grid.inside       = 1:size(cfg.grid.pos, 1);
        cfg.grid.outside      = [];
        cfg.lcmv.keepfilter   = 'yes';
        source_idx            = ft_sourceanalysis(cfg, timelock);
        beamformer            = source_idx.avg.filter{1};
        
        source_timecourses.filter = [source_timecourses.filter; {beamformer}];
        
        % compute mom (source power in each direction)
        pow_data = [];
        pow_data.label = {'pow_x', 'pow_y', 'pow_z'};
        pow_data.time = data_continuous_concat_down.time;
        for i=1:length(data_continuous_concat_down.trial) 
            pow_data.trial{i} = beamformer * data_continuous_concat_down.trial{i}(chansel,:);
        end
        
        source_timecourses.mom = [source_timecourses.mom; {pow_data.trial{i}}];

        % project the source to its strongest orientation, i.e. the 
        % direction that explains most of the source variance. This is done
        % by single-value decomposition that is equivalent to PCA
        mom_voxel = cat(2, pow_data.trial{:});
        [u, s, v] = svd(mom_voxel, 'econ');

        % this is equal to the first column of matrix V, apart from the scaling with s(1,1)
        timeseries_voxel = u(:,1)' * mom_voxel;
        source_timecourses.trial = [source_timecourses.trial; timeseries_voxel];
    
    else
        source_timecourses.trial = [source_timecourses.trial; zeros(1,length(data_continuous_concat_down.time{1}))];
    end
    
end

source_timecourses.trial = {source_timecourses.trial};

%% t = toc;
%% fprintf(['\n###############\nFinished - calc_source_timecourses took ' num2str(t) ' seconds.\n\n'])

end

