%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% BATCHING Power Spectrum Analysis
%%% 
%%%  1) Compute power analysis for each subject and each run separatly
%%%     ouput: one .mat per run / 102channels * 40freq
%%%
%%%  2) Hight_contact : Compute power analysis for each condition separatly (average the two runs over all the subjects)
%%%     ouput: one .mat / 102channels * 40freq
%%%
%%%  3) Low_contact : Compute power analysis for each condition separatly (average the two runs over all the subjects)
%%%     ouput: one .mat / 102channels * 40freq
%%%
%%%  4) Compute power Spectrum for each suject and each run taken separatly
%%%     ouput: 4 .mat (1 per condition and run) / 20sujets * 102channels * 40freq
%%%
%%%  5) Compute power Spectrum for each suject and each condition taken separatly (average the runs)
%%%     ouput: 4 .mat (1 per condition and run) / 20sujets * 102channels * 40freq
%%% 
%%%  on magnometers channels
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Compute power analysis for each subject and each run separatly
%ouput: one .mat per run / 102channels * 40freq

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
            Power_Analysis(input_filename, data_no_muscle, condition) %run power analysis on clean Mag sensors
            end
	    end

  catch
	disp(['Something was wrong with Subject' int2str(i) 'for Power Analysis! Continuing with next in line']);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Hight_contact : Compute power analysis for each condition separatly (average the two runs over all the subjects)
%ouput: one .mat / 102channels * 40freq

%Load and average the data (high_contact)
clear all
subjects = [1, 2, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 19, 20, 21, 22];  

name = genvarname(repmat({'data'},38,1))
n=1

for i = subjects
  try
        eval(['subject' int2str(i)])
        for S = [1 2]  %Pour les 2 sessions
	    datadir =['subjectdata.datadir' int2str(S)]
        	for R = [1 2]  %Pour les 2 runs
                if i == 1 || i == 6 || i == 10 || i == 11 || i == 13 || i == 15 || i == 16 || i == 20 || i == 22;
                input_filename = [subjectdata.subjectdir filesep eval(datadir) filesep 'run0' int2str(R) '_tsss_realign_ica_eye_corrected_grad_ica_cardio_corrected_grad_power.mat']
                else input_filename = [subjectdata.subjectdir filesep eval(datadir) filesep 'run0' int2str(R) '_tsss_ica_eye_corrected_grad_ica_cardio_corrected_grad_power.mat']
                end 
            condition = eval(['subjectdata.condition' int2str(S) int2str(R)])
            if strcmp(condition,'High_contact')
            load(input_filename)
            name{n} = power
            n = n+1
            end
            end
	    end

  catch
	disp(['Something was wrong with Subject' int2str(i) 'for Power Analysis! Continuing with next in line']);
  end
end

%%Average the data

cfg = []
%cfg.keepindividual = 'yes' or 'no' (default = 'no')
%cfg.foilim         = [fmin fmax] or 'all', to specify a subset of frequencies (default = 'all')
%cfg.toilim         = [tmin tmax] or 'all', to specify a subset of latencies (default = 'all')
cfg.channel = 'MEG*1';
%cfg.parameter      = string or cell-array of strings indicating which parameter(s) to average. default is set to 'powspctrm', if it is present in the data.
power_high_contact = ft_freqgrandaverage(cfg, name{:})

%Save the output 
filename = [ '/pclxserver/raid11/data/SELFI/Grand_average/High_contact_power.mat' ] ;
display( [ 'Saving ' filename ] )
save( filename, 'power_high_contact' );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Low_contact : Compute power analysis for each condition separatly (average the two runs over all the subjects)
%ouput: one .mat / 102channels * 40freq

%Load and average the data (Low_contact)

clear all
subjects = [1, 2, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 19, 20, 21, 22];  
name = genvarname(repmat({'data'},38,1))
n=1

for i = subjects
  try
        eval(['subject' int2str(i)])
        for S = [1 2]  %Pour les 2 sessions
	    datadir =['subjectdata.datadir' int2str(S)]
        	for R = [1 2]  %Pour les 2 runs
                if i == 1 || i == 6 || i == 10 || i == 11 || i == 13 || i == 15 || i == 16 || i == 20 || i == 22;
                input_filename = [subjectdata.subjectdir filesep eval(datadir) filesep 'run0' int2str(R) '_tsss_realign_ica_eye_corrected_grad_ica_cardio_corrected_grad_power.mat']
                else input_filename = [subjectdata.subjectdir filesep eval(datadir) filesep 'run0' int2str(R) '_tsss_ica_eye_corrected_grad_ica_cardio_corrected_grad_power.mat']
                end 
            condition = eval(['subjectdata.condition' int2str(S) int2str(R)])
            if strcmp(condition,'Low_contact')
            load(input_filename)
            name{n} = power
            n = n+1
            end
            end
	    end

  catch
	disp(['Something was wrong with Subject' int2str(i) 'for Power Analysis! Continuing with next in line']);
  end
end

%%Average the data

cfg = []
%cfg.keepindividual = 'yes' or 'no' (default = 'no')
%cfg.foilim         = [fmin fmax] or 'all', to specify a subset of frequencies (default = 'all')
%cfg.toilim         = [tmin tmax] or 'all', to specify a subset of latencies (default = 'all')
cfg.channel = 'MEG*1';
%cfg.parameter      = string or cell-array of strings indicating which parameter(s) to average. default is set to 'powspctrm', if it is present in the data.
power_low_contact = ft_freqgrandaverage(cfg, name{:})

%Save the output 
filename = [ '/pclxserver/raid11/data/SELFI/Grand_average/Low_contact_power.mat' ] ;
display( [ 'Saving ' filename ] )
save( filename, 'power_low_contact' );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Compute power Spectrum for each suject and each run taken separatly
%ouput: 4 .mat (1 per condition and run) / 20sujets * 102channels * 40freq

clear all
subjects = [1, 2, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 19, 20, 21, 22];  

allsubj_power_High_contact_R1 = genvarname(repmat({'data'},19,1))
a=1

allsubj_power_High_contact_R2 = genvarname(repmat({'data'},19,1))
b=1

allsubj_power_Low_contact_R1 = genvarname(repmat({'data'},19,1))
c=1

allsubj_power_Low_contact_R2 = genvarname(repmat({'data'},19,1))
d=1


for i = subjects
  try
        eval(['subject' int2str(i)])
        for S = [1 2]  %Pour les 2 sessions
	    datadir =['subjectdata.datadir' int2str(S)]
        	for R = [1 2]  %Pour les 2 runs
                if i == 1 || i == 6 || i == 10 || i == 11 || i == 13 || i == 15 || i == 16 || i == 20 || i == 22;
                input_filename = [subjectdata.subjectdir filesep eval(datadir) filesep 'run0' int2str(R) '_tsss_realign_ica_eye_corrected_grad_ica_cardio_corrected_grad_power.mat']
                else input_filename = [subjectdata.subjectdir filesep eval(datadir) filesep 'run0' int2str(R) '_tsss_ica_eye_corrected_grad_ica_cardio_corrected_grad_power.mat']
                end 
            condition = eval(['subjectdata.condition' int2str(S) int2str(R)])
            load(input_filename)
            
            if strcmp(condition,'High_contact') & R == 1
            allsubj_power_High_contact_R1{a} = power
            a = a+1    
            elseif  strcmp(condition,'High_contact') & R == 2
            allsubj_power_High_contact_R2{b} = power
            b = b+1   
            elseif  strcmp(condition,'Low_contact') & R == 1
            allsubj_power_Low_contact_R1{c} = power
            c = c+1 
            elseif  strcmp(condition,'Low_contact') & R == 2
            allsubj_power_Low_contact_R2{d} = power
            d = d+1 
            end
            end
	    end

  catch
	disp(['Something was wrong with Subject' int2str(i) 'for Power Analysis! Continuing with next in line']);
  end
end

% Save the output 
save( '/pclxserver/raid11/data/SELFI/Grand_average/allsubj_power_High_contact_R1.mat', 'allsubj_power_High_contact_R1' );
save( '/pclxserver/raid11/data/SELFI/Grand_average/allsubj_power_High_contact_R2.mat', 'allsubj_power_High_contact_R2' );
save( '/pclxserver/raid11/data/SELFI/Grand_average/allsubj_power_Low_contact_R1.mat', 'allsubj_power_Low_contact_R1' );
save( '/pclxserver/raid11/data/SELFI/Grand_average/allsubj_power_Low_contact_R2.mat', 'allsubj_power_Low_contact_R2' );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Compute power Spectrum for each suject and each condition taken separatly (average the runs)
%ouput: 4 .mat (1 per condition and run) / 20sujets * 102channels * 40freq

clear all
subjects = [1, 2, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 19, 20, 21, 22];  

allsubj_power_High_contact = genvarname(repmat({'data'},19,1))
a=1

allsubj_power_Low_contact = genvarname(repmat({'data'},19,1))
c=1


for i = subjects
  try
        eval(['subject' int2str(i)])
        for S = [1 2]  %Pour les 2 sessions
	    datadir =['subjectdata.datadir' int2str(S)]
        	for R = [1 2]  %Pour les 2 runs
                if i == 1 || i == 6 || i == 10 || i == 11 || i == 13 || i == 15 || i == 16 || i == 20 || i == 22;
                input_filename = [subjectdata.subjectdir filesep eval(datadir) filesep 'run0' int2str(R) '_tsss_realign_ica_eye_corrected_grad_ica_cardio_corrected_grad_power.mat']
                else input_filename = [subjectdata.subjectdir filesep eval(datadir) filesep 'run0' int2str(R) '_tsss_ica_eye_corrected_grad_ica_cardio_corrected_grad_power.mat']
                end 
            condition = eval(['subjectdata.condition' int2str(S) int2str(R)])
            load(input_filename)
            
            if strcmp(condition,'High_contact') 
            allsubj_power_High_contact{a} = power
            a = a+1       
            elseif  strcmp(condition,'Low_contact') 
            allsubj_power_Low_contact{c} = power
            c = c+1  
            end
            end
	    end

  catch
	disp(['Something was wrong with Subject' int2str(i) 'for Power Analysis! Continuing with next in line']);
  end
end

% Save the output 
save( '/pclxserver/raid11/data/SELFI/Grand_average/allsubj_power_High_contact.mat', 'allsubj_power_High_contact' );
save( '/pclxserver/raid11/data/SELFI/Grand_average/allsubj_power_Low_contact.mat', 'allsubj_power_Low_contact' );






