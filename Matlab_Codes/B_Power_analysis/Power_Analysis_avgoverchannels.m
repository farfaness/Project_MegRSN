function [] = Power_Analysis_avgoverchannels(input_filename, data_no_muscle, condition)
           
%% Define trials only on clean data segments

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

%% Compute Power
cfg = [];
cfg.channel = 'MEG*1';
cfg.keeptrials = 'no' %return individual trials or average ('no')
cfg.method     = 'mtmfft' %analyses an entire spectrum for the entire data length, implements multitaper frequency transformation %or 'mtmconvol' 'wavelet' 'tfr' 'mvar'
cfg.taper      = 'hanning' 
cfg.foi     = [1:1:40];   %freq band of interest with step of 1Hz (2nd nbr)
cfg.output     = 'pow' % 'pow' return the power-spectra,'powandcsd' return the power and the cross-spectra, 'fourier'   return the complex Fourier-spectra
%cfg.pad        = %number, 'nextpow2', or 'maxperlen' (default), length in seconds to which the data can be padded out. Padding will determine your spectral resolution. 
%cfg.padtype     = %string, type of padding (default 'zero', see ft_preproc_padding)
cfg.tapsmofrq  = 0.5 %number, the amount of spectral smoothing through multi-tapering. Note that 4 Hz smoothing means plus-minus 4 Hz, i.e. a 8 Hz smoothing box.
power = ft_freqanalysis(cfg, data_trial)

%% Average over channels
cfg = [];
cfg.channel = 'MEG*1';
cfg.avgoverchan = 'yes';
% cfg.frequency   = [beg end]
% cfg.avgoverfreq = string, can be 'yes' or 'no' (default = 'no')
power_avgoverchannels = ft_selectdata(cfg, power)


%% Save the output 
[path,name,] = fileparts(input_filename) ;
Power_filename = [ path '/' name '_power_avgoverchannels.mat' ] ;
display( [ 'Saving ' Power_filename ] )
save( Power_filename, 'power_avgoverchannels');


end