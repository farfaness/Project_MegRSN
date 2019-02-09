%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SCRIPT_Visu
%%% 
%%% Script to visualise the data : 
%%%   1) before ica correction
%%%   2) after ica blink correction 
%%%   3) after ica cardio correction
%%% 
%%% on magnometers, and averaged trials of 8s 
%%%
%%% Useful to compare the data before and after ica correction 
%%% to check if the correction was done correctly
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Start by defining the inputs:

% clear all; subject_number= 2 ; session_number= 1; run_number = 1
function [] = visualize_before_and_after_ica_correction(subject_number, session_number, run_number)
eval(['subject' int2str(subject_number)])

if isequal(session_number, 1)
   datadir = subjectdata.datadir1
elseif isequal(session_number, 2)
   datadir = subjectdata.datadir2
else
   disp('Wrong session number')
end

if isequal(run_number, 1)
   run = 'run01'
elseif isequal(run_number, 2)
   run = 'run02'
else
   disp('Wrong run number')
end

input_filename_tsss = [subjectdata.subjectdir filesep datadir filesep run '_tsss.fif']
input_filename_ica_eye = [subjectdata.subjectdir filesep datadir filesep run '_tsss_ica_eye_corrected.fif']
input_filename_ica_cardio = [subjectdata.subjectdir filesep datadir filesep run '_tsss_ica_eye_corrected_ica_cardio_corrected.fif']

data_no_muscle = eval(['subjectdata.data_no_muscle' int2str(session_number) int2str(run_number)])
A= size(data_no_muscle)
nb_artefact = A(1)-1

artefact= NaN(nb_artefact,2)

for i = 1:nb_artefact
    artefact (i, 1) = data_no_muscle (i, 2)
    artefact (i, 2) = data_no_muscle ((i+1), 1)
    i = i+1
end

%% Preprocess the data before and after the ICA correction for visual inspection (define and average 37 trial of 8s = 296s)

%Data before correction
cfg         = [];
cfg.dataset = input_filename_tsss ;
cfg.trialfun             = 'ft_trialfun_general';
cfg.trialdef.triallength = 8;                   % in seconds
cfg.trialdef.ntrials     = 37;                 % i.e. the complete file : 37 * 8s = 296s
cfg = ft_definetrial(cfg);    %The trial definition "trl" (Nx3 matrix, N: nbr of trials, offset of 0: the first sample of the trial corresponds to the trigger) is added to the output configuration structure

cfg.artfctdef.muscle.artifact = artefact
cfg.artfctdef.reject          = 'complete'
cfg.artfctdef.feedback        = 'yes' %For a visual feedback on the artefacts
cfg = ft_rejectartifact(cfg);

cfg.channel = 'MEG*1';
data_tsss = ft_preprocessing(cfg); %The info in "trl" is put into 2 differents field of the datastructure: data.sampleinfo (start end) and data.trialinfo (stimulus code)

%Data after ica blink correction
cfg         = [];
cfg.dataset = input_filename_ica_eye ;
cfg.trialfun             = 'ft_trialfun_general';
cfg.trialdef.triallength = 8;                   % in seconds
cfg.trialdef.ntrials     = 37;                 % i.e. the complete file : 37 * 8s = 296s
cfg = ft_definetrial(cfg);    %The trial definition "trl" (Nx3 matrix, N: nbr of trials, offset of 0: the first sample of the trial corresponds to the trigger) is added to the output configuration structure

cfg.artfctdef.muscle.artifact = artefact
cfg.artfctdef.reject          = 'complete'
cfg = ft_rejectartifact(cfg);

cfg.channel = 'MEG*1';
data_ica_eye = ft_preprocessing(cfg); %The info in "trl" is put into 2 differents field of the datastructure: data.sampleinfo (start end) and data.trialinfo (stimulus code)


%Data after ica cardio correction
cfg         = [];
cfg.dataset = input_filename_ica_cardio ;
cfg.trialfun             = 'ft_trialfun_general';
cfg.trialdef.triallength = 8;                   % in seconds
cfg.trialdef.ntrials     = 37;                 % i.e. the complete file : 37 * 8s = 296s
cfg = ft_definetrial(cfg);    %The trial definition "trl" (Nx3 matrix, N: nbr of trials, offset of 0: the first sample of the trial corresponds to the trigger) is added to the output configuration structure

cfg.artfctdef.muscle.artifact = artefact
cfg.artfctdef.reject          = 'complete'
cfg = ft_rejectartifact(cfg);

cfg.channel = 'MEG*1';
data_ica_cardio = ft_preprocessing(cfg); %The info in "trl" is put into 2 differents field of the datastructure: data.sampleinfo (start end) and data.trialinfo (stimulus code)


%% plot the data before and after the ICA correction for visual inspection

cfg = [];
cfg.xlim = [0 8];
cfg.ylim = [-1e-13 3e-13];
cfg.channel = 'MEG*1';
figure; ft_singleplotER(cfg, data_tsss, data_ica_eye, data_ica_cardio);
end


