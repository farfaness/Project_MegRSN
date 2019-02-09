%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% BATCHING
%%% 
%%% ICA to correct for blinks first,
%%%   and then another ICA to correct for cardio artefacts separatly 
%%%   on the file corrected for blinks
%%% 
%%% on magnometers channels
%%%
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% ICA to correct blinks first 

clear all
number_of_subjects = 22  

%Blinks correction

for i = 1:number_of_subjects
   try
        eval(['subject' int2str(i)])
        for S = [1]  %Pour les 2 sessions
	    datadir =['subjectdata.datadir' int2str(S)]
        	for R = [1]  %Pour les 2 runs
        	input_filename = [subjectdata.subjectdir filesep eval(datadir) filesep 'run0' int2str(R) '_tsss.fif']
            data_no_muscle = eval(['subjectdata.data_no_muscle' int2str(S) int2str(R)])
            do_ica_blinks_mag(input_filename, data_no_muscle) %run ICA on MEG data
        	find_all_components_blinks_mag(input_filename, data_no_muscle)  %correlation between bio 2 and the components
         	[path,name,~] = fileparts(input_filename) ;
        	components_filename = [ path '/' name '_ica_comp_blinks.mat' ];
        	components_to_remove = [ path '/' name '_components_to_remove_blinks.mat' ]
        	reject_all_components_blinks_mag(input_filename, components_filename, components_to_remove ) %rejects the components that correlate
        	end
	end

  catch
 	disp(['Something was wrong with Subject' int2str(i) 'for cardio ICA! Continuing with next in line']);
     
   end
end

%%%%% ICA to correct for Cardio artefact on the fif file corrected for blinks 

clear all
number_of_subjects = 22  

for i = 1:number_of_subjects
  try
        eval(['subject' int2str(i)])
        for S = [1 2]  %Pour les 2 sessions
	    datadir =['subjectdata.datadir' int2str(S)]
        	for R = [1 2]  %Pour les 2 runs
        	input_filename = [subjectdata.subjectdir filesep eval(datadir) filesep 'run0' int2str(R) '_tsss_ica_eye_corrected.fif']
            data_no_muscle = eval(['subjectdata.data_no_muscle' int2str(S) int2str(R)])
        	do_ica_cardio_mag(input_filename, data_no_muscle) %run ICA on MEG data
        	find_all_components_cardio_mag(input_filename, data_no_muscle)  %correlation between bio 3 and the components
        	[path,name,~] = fileparts(input_filename) ;
        	components_filename = [ path '/' name '_ica_comp_cardio.mat' ];
        	components_to_remove = [ path '/' name '_components_to_remove_cardio.mat' ]
        	reject_all_components_cardio_mag(input_filename, components_filename, components_to_remove ) %rejects the components that correlate
        	end
	end

  catch
	disp(['Something was wrong with Subject' int2str(i) 'for cardio ICA! Continuing with next in line']);
    
  end
end
