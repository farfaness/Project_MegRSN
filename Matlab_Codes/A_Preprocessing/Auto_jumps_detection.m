%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SCRIPT_Auto_jumps_detection
%%% 
%%% Script to detect jumps in the data automaticaly 
%%% 
%%% on all the channels and all the trials
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% load the data
cfg         = [];
cfg.dataset = input_filename ;
cfg.channel = 'MEG';
cfg.trialfun             = 'ft_trialfun_general';
cfg.trialdef.triallength = 1;                   % in seconds
cfg.trialdef.ntrials     = inf;                 % i.e. the complete file
cfg = ft_definetrial(cfg);                      % this creates 1-second data segments
trl                    = cfg.trl(1:301,:); % we'll use our 301 trials

data        = ft_preprocessing(cfg);

% jump
cfg                    = [];
cfg.trl = trl;
%cfg.datafile   = input_filename; % pas necessaire lorsque data en input de ft_artifact_zvalue
%cfg.headerfile = input_filename; % pas necessaire lorsque data en input de ft_artifact_zvalue
cfg.continuous = 'yes';
 
% channel selection, cutoff and padding
cfg.artfctdef.zvalue.channel    = 'MEG';
cfg.artfctdef.zvalue.cutoff     = 20;
cfg.artfctdef.zvalue.trlpadding = 0;
cfg.artfctdef.zvalue.artpadding = 0;
cfg.artfctdef.zvalue.fltpadding = 0;
 
% algorithmic parameters
cfg.artfctdef.zvalue.cumulative    = 'yes';
cfg.artfctdef.zvalue.medianfilter  = 'yes';
cfg.artfctdef.zvalue.medianfiltord = 9;
cfg.artfctdef.zvalue.absdiff       = 'yes';
 
% make the process interactive
cfg.artfctdef.zvalue.interactive = 'yes';
 
[cfg, artifact_jump] = ft_artifact_zvalue(cfg, data);

