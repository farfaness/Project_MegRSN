%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SCRIPT_Batch_BAND_PASS_filtering_epoched
%%%
%%% This script filter the data in alpha, beta and theta frequency bands
%%% for the two conditions dataset, on epoched data
%%%
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% ALpha Band-pass filtering off the epoched data (8 seconds clean data trials)

clear all
close all
 
subjects = [1, 2, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 19, 20, 21, 22]; 
lpfreq = 8
hpfreq = 13
name_band_freq = 'alpha'
 
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
             Band_pass_filtering_epoched_data(input_filename, data_no_muscle, condition, lpfreq, hpfreq, name_band_freq) 
             end
 	    end
 
   catch
 	disp(['Something was wrong with Subject' int2str(i) 'for Epoched Frequency Filtering! Continuing with next in line']);
   end
 end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%% Beta Band-pass filtering off the epoched data (8 seconds clean data trials)

clear all
close all
 
subjects = [1, 2, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 19, 20, 21, 22]; 
lpfreq = 13
hpfreq = 30
name_band_freq = 'beta'
 
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
             Band_pass_filtering_epoched_data(input_filename, data_no_muscle, condition, lpfreq, hpfreq, name_band_freq) 
             end
 	    end
 
   catch
 	disp(['Something was wrong with Subject' int2str(i) 'for Epoched Frequency Filtering! Continuing with next in line']);
   end
 end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Theta Band-pass filtering off the epoched data (8 seconds clean data trials)

clear all
close all

subjects = [1, 2, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 19, 20, 21, 22]; 
lpfreq = 4
hpfreq = 7
name_band_freq = 'theta'

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
            Band_pass_filtering_epoched_data(input_filename, data_no_muscle, condition, lpfreq, hpfreq, name_band_freq) 
            end
	    end

  catch
	disp(['Something was wrong with Subject' int2str(i) 'for Epoched Frequency Filtering! Continuing with next in line']);
  end
end

