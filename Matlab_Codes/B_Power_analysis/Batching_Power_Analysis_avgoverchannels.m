%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% BATCHING Power Spectrum Analysis Averaged over channels
%%% 
%%%  1) Compute power analysis for each subject and each run separatly on averaged channels
%%%     ouput: one .mat per run / 1 channel * 40freq
%%%
%%%  2) Hight_contact : Compute power analysis for each condition separatly (the two runs over all the subjects in one file)
%%%     ouput: one .mat / 1 channel * 40freq
%%%
%%%  3) Low_contact : Compute power analysis for each condition separatly (the two runs over all the subjects in one file)
%%%     ouput: one .mat / 1 channel * 40freq
%%%
%%%  4) Compute power Spectrum for each suject and each run taken separatly
%%%     ouput: 4 .mat (1 per condition and run) / 20sujets * 1 channel * 40freq
%%%
%%%  5) Compute power Spectrum for each suject and each condition taken separatly (average the runs)
%%%     ouput: 4 .mat (1 per condition and run) / 20sujets * 1 channel * 40freq
%%% 
%%%  on magnometers channels
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Compute power analysis for each subject and for each run separatly on averaged channels
%ouput: one .mat per run / 1channel * 40freq

clear all
subjects = [1, 2, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 19, 20, 21, 22];  

for i = subjects
  try
        eval(['subject' int2str(i)])
        for S = [1 2]  %Pour les 2 sessions
	    datadir =['subjectdata.datadir' int2str(S)]
        	for R = [1 2]  %Pour les 2 runs
                if i == 1 || i == 6 || i == 10 || i == 11 || i == 13 || i == 15 || i == 16 || i == 20 || i == 22;
                input_filename = [subjectdata.subjectdir filesep eval(datadir) filesep 'run0' int2str(R) '_tsss_realign_ica_eye_corrected_grad_ica_cardio_corrected_grad.fif']
                else input_filename = [subjectdata.subjectdir filesep eval(datadir) filesep 'run0' int2str(R) '_tsss_ica_eye_corrected_grad_ica_cardio_corrected_grad.fif']
                end 
            data_no_muscle = eval(['subjectdata.data_no_muscle' int2str(S) int2str(R)])
            condition = eval(['subjectdata.condition' int2str(S) int2str(R)])
            Power_Analysis_avgoverchannels(input_filename, data_no_muscle, condition)
            end
	    end

  catch
	disp(['Something was wrong with Subject' int2str(i) 'for Power Analysis! Continuing with next in line']);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Hight_contact : Compute power analysis for each condition separatly (the two runs over all the subjects in one file)
%ouput: one .mat / 1 channel * 40freq

%Load the data
clear all
input_filename = '/pclxserver/raid11/data/SELFI/Grand_average/High_contact_power.mat'
load (input_filename);

%Average the channels
cfg = [];
cfg.channel = 'MEG*1';
cfg.avgoverchan = 'yes';
% cfg.frequency   = [beg end]
% cfg.avgoverfreq = string, can be 'yes' or 'no' (default = 'no')
power_high_contact_avgoverchannels = ft_selectdata(cfg, power_high_contact)

%Save the output 
filename = [ '/pclxserver/raid11/data/SELFI/Grand_average/High_contact_power_avgoverchannels.mat' ] ;
display( [ 'Saving ' filename ] )
save( filename, 'power_high_contact_avgoverchannels' );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Low_contact : Compute power analysis for each condition separatly (the two runs over all the subjects in one file)
%ouput: one .mat / 1 channel * 40freq

%Load and average the data (Low_contact)

clear all

%Load the data
clear all
input_filename = '/pclxserver/raid11/data/SELFI/Grand_average/Low_contact_power.mat'
load (input_filename);

%Average the channels
cfg = [];
cfg.channel = 'MEG*1';
cfg.avgoverchan = 'yes';
% cfg.frequency   = [beg end]
% cfg.avgoverfreq = string, can be 'yes' or 'no' (default = 'no')
power_low_contact_avgoverchannels = ft_selectdata(cfg, power_low_contact)

%Save the output 
filename = [ '/pclxserver/raid11/data/SELFI/Grand_average/Low_contact_power_avgoverchannels.mat' ] ;
display( [ 'Saving ' filename ] )
save( filename, 'power_low_contact_avgoverchannels' );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Compute power Spectrum for each suject and each run taken separatly
%ouput: 4 .mat (1 per condition and run) / 20sujets * 1 channel * 40freq

clear all
subjects = [1, 2, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 19, 20, 21, 22];

allsubj_power_High_contact_R1_avgoverchannels = genvarname(repmat({'data'},19,1))
a=1

allsubj_power_High_contact_R2_avgoverchannels = genvarname(repmat({'data'},19,1))
b=1

allsubj_power_Low_contact_R1_avgoverchannels = genvarname(repmat({'data'},19,1))
c=1

allsubj_power_Low_contact_R2_avgoverchannels = genvarname(repmat({'data'},19,1))
d=1


for i = subjects
  try
        eval(['subject' int2str(i)])
        for S = [1 2]  %Pour les 2 sessions
	    datadir =['subjectdata.datadir' int2str(S)]
        	for R = [1 2]  %Pour les 2 runs
                if i == 1 || i == 6 || i == 10 || i == 11 || i == 13 || i == 15 || i == 16 || i == 20 || i == 22;
                input_filename = [subjectdata.subjectdir filesep eval(datadir) filesep 'run0' int2str(R) '_tsss_realign_ica_eye_corrected_grad_ica_cardio_corrected_grad_power_avgoverchannels.mat']
                else input_filename = [subjectdata.subjectdir filesep eval(datadir) filesep 'run0' int2str(R) '_tsss_ica_eye_corrected_grad_ica_cardio_corrected_grad_power_avgoverchannels.mat']
                end 
            condition = eval(['subjectdata.condition' int2str(S) int2str(R)])
            load(input_filename)
            
            if strcmp(condition,'High_contact') & R == 1
            allsubj_power_High_contact_R1_avgoverchannels{a} = power_avgoverchannels
            a = a+1    
            elseif  strcmp(condition,'High_contact') & R == 2
            allsubj_power_High_contact_R2_avgoverchannels{b} = power_avgoverchannels
            b = b+1   
            elseif  strcmp(condition,'Low_contact') & R == 1
            allsubj_power_Low_contact_R1_avgoverchannels{c} = power_avgoverchannels
            c = c+1 
            elseif  strcmp(condition,'Low_contact') & R == 2
            allsubj_power_Low_contact_R2_avgoverchannels{d} = power_avgoverchannels
            d = d+1 
            end
            end
	    end

  catch
	disp(['Something was wrong with Subject' int2str(i) 'for Power Analysis! Continuing with next in line']);
  end
end

% Save the output 
save( '/pclxserver/raid11/data/SELFI/Grand_average/allsubj_power_High_contact_R1_avgoverchannels.mat', 'allsubj_power_High_contact_R1_avgoverchannels' );
save( '/pclxserver/raid11/data/SELFI/Grand_average/allsubj_power_High_contact_R2_avgoverchannels.mat', 'allsubj_power_High_contact_R2_avgoverchannels' );
save( '/pclxserver/raid11/data/SELFI/Grand_average/allsubj_power_Low_contact_R1_avgoverchannels.mat', 'allsubj_power_Low_contact_R1_avgoverchannels' );
save( '/pclxserver/raid11/data/SELFI/Grand_average/allsubj_power_Low_contact_R2_avgoverchannels.mat', 'allsubj_power_Low_contact_R2_avgoverchannels' );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Compute power Spectrum for each suject and each condition taken separatly (average the runs)
%ouput: 4 .mat (1 per condition and run) / 20sujets * 1 channel * 40freq

clear all
subjects = [1, 2, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 19, 20, 21, 22];

allsubj_power_High_contact_avgoverchannels = genvarname(repmat({'data'},38,1))
a=1

allsubj_power_Low_contact_avgoverchannels = genvarname(repmat({'data'},38,1))
c=1


for i = subjects
  try
        eval(['subject' int2str(i)])
        for S = [1 2]  %Pour les 2 sessions
	    datadir =['subjectdata.datadir' int2str(S)]
        	for R = [1 2]  %Pour les 2 runs
                if i == 1 || i == 6 || i == 10 || i == 11 || i == 13 || i == 15 || i == 16 || i == 20 || i == 22;
                input_filename = [subjectdata.subjectdir filesep eval(datadir) filesep 'run0' int2str(R) '_tsss_realign_ica_eye_corrected_grad_ica_cardio_corrected_grad_power_avgoverchannels.mat']
                else input_filename = [subjectdata.subjectdir filesep eval(datadir) filesep 'run0' int2str(R) '_tsss_ica_eye_corrected_grad_ica_cardio_corrected_grad_power_avgoverchannels.mat']
                end 
            condition = eval(['subjectdata.condition' int2str(S) int2str(R)])
            load(input_filename)
            
            if strcmp(condition,'High_contact') 
            allsubj_power_High_contact_avgoverchannels{a} = power_avgoverchannels
            a = a+1       
            elseif  strcmp(condition,'Low_contact') 
            allsubj_power_Low_contact_avgoverchannels{c} = power_avgoverchannels
            c = c+1  
            end
            end
	    end

  catch
	disp(['Something was wrong with Subject' int2str(i) 'for Power Analysis! Continuing with next in line']);
  end
end

% Average the two runs
y = 1
x= 1
for y = 1:19
x = y * 2    
cfg = []
cfg.parameter  = 'powspctrm'
allsubj_power_High_contact_avgoverchannels{y, 1} = ft_freqgrandaverage(cfg, allsubj_power_High_contact_avgoverchannels{x-1}, allsubj_power_High_contact_avgoverchannels{x})
y = y + 1
end

z = 1
v = 1
for z = 1:19
v = z * 2
cfg = []
cfg.parameter  = 'powspctrm'
allsubj_power_Low_contact_avgoverchannels{z, 1} = ft_freqgrandaverage(cfg, allsubj_power_Low_contact_avgoverchannels{v-1}, allsubj_power_Low_contact_avgoverchannels{v})
z = z+1
end

% Save the output 
save( '/pclxserver/raid11/data/SELFI/Grand_average/allsubj_power_High_contact_avgoverchannels.mat', 'allsubj_power_High_contact_avgoverchannels' );
save( '/pclxserver/raid11/data/SELFI/Grand_average/allsubj_power_Low_contact_avgoverchannels.mat', 'allsubj_power_Low_contact_avgoverchannels' );

