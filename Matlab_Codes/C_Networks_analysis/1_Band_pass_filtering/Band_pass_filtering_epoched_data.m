function [] = Band_pass_filtering_epoched_data( input_filename, data_no_muscle, condition, lpfreq, hpfreq, name_band_freq )

%%%%% BAND_PASS_FILTERING 

%   Band_pass filter data in fonction of the specified band of frequency
%   Define trials only on clean data segments

%Determine condition : 1 for High Contact, 2 for Low Contact
if strcmp(condition,'High_contact')
     cond = 1;
elseif strcmp(condition,'Low_contact')
     cond = 2;
end

%Define the trials
A= size(data_no_muscle)
nb_artefact = A(1)-1
artefact= NaN(nb_artefact,2)

for i = 1:nb_artefact
    artefact (i, 1) = data_no_muscle (i, 2)
    artefact (i, 2) = data_no_muscle ((i+1), 1)
    i = i+1
end

cfg         = [];
cfg.dataset = input_filename ;
cfg.trialfun             = 'ft_trialfun_general';
cfg.trialdef.triallength = 8;                   % in seconds
cfg.trialdef.ntrials     = 37;                 % i.e. the complete file
cfg = ft_definetrial(cfg);    %The trial definition "trl" (Nx3 matrix, N: nbr of trials, offset of 0: the first sample of the trial corresponds to the trigger) is added to the output configuration structure

cfg.artfctdef.muscle.artifact = artefact
cfg.artfctdef.reject          = 'complete'
%cfg.artfctdef.feedback        = 'yes' %For a visual feedback of the trials rejected
cfg = ft_rejectartifact(cfg);

data_trial = ft_preprocessing(cfg); %The info in "trl" is put into 2 differents field of the datastructure: data.sampleinfo (start end) and data.trialinfo (stimulus code)

data_trial.trialinfo = cond

%% compute Band-pass filtering the epoched data

cfg = [];
cfg.lpfilter      = 'yes';
cfg.hpfilter      = 'yes';
cfg.lpfreq        = lpfreq;
cfg.hpfreq        = hpfreq;
%cfg.demean        = 'yes';  % To do for the covariance matrix (not necessary because ft_timelock automaticaly do it)
[data_epoched_filtered] = ft_preprocessing(cfg, data_trial);

%% Save the output 
[path,name,] = fileparts(input_filename) ;
Band_pass_filtered_filename = [ path '/' name '_epoched_' name_band_freq '.mat' ] ;
display( [ 'Saving ' Band_pass_filtered_filename  ] )
save( Band_pass_filtered_filename, 'data_epoched_filtered' );


end

