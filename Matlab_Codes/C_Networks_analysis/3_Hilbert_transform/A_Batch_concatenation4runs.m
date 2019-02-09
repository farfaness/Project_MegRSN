%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Batch_concatenation4runs
%%%
%%%  For each subject : 
%%%  - takes the source time courses
%%%  - hilbert transform
%%%  - downsample
%%%  - reject bad segments
%%%  - concatenates the 4 runs for each subject
%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

% Indicate the Frequency band of interest:
name_band_freq = 'theta';

subjects = [1, 2, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 19, 20, 21, 22];    

 for i = subjects
     
   eval(['subject' int2str(i)])
   Outputdir = eval('subjectdata.subjectdir');
   
   % Load the source_timecourse
   disp('Loading data...');
   load([subjectdata.subjectdir filesep 'HC_and_LC_' name_band_freq '_source_timecourses1cm'], ['source_timecourses1cm'])
       
   % Hilbert transform
   alpha_source_timecourses_hilbert =  source_timecourses1cm;
   alpha_source_timecourses_hilbert.trial{1} = ft_preproc_hilbert(source_timecourses1cm.trial{1}, 'abs');
   
   %%% Not saved for disk space reasons
   %save([subjectdata.subjectdir filesep 'HC_and_LC_' name_band_freq '_source_timecourses_hilbert'], 'alpha_source_timecourses_hilbert', '-v7.3')
   %load([subjectdata.subjectdir filesep 'HC_and_LC_' name_band_freq '_source_timecourses_hilbert'])
 
   VERIF_time = size(alpha_source_timecourses_hilbert.trial{1, 1});
 
   %% Downsampling
   disp('Downsampling...')
   % "dummy" channel labels for each grid have to be created for
   % ft_resampledata to accept the input data
   alpha_source_timecourses_hilbert_dummy_bis = alpha_source_timecourses_hilbert;  
   alpha_source_timecourses_hilbert_dummy_bis.label = create_random_labels(size(alpha_source_timecourses_hilbert.trial{1},1));
   
   cfg = [];
   cfg.resamplefs = 1;     
   cfg.detrend    = 'no';
   cfg.demean     = 'no';
   alpha_source_timecourses_hilbert_downsampled_bis = ft_resampledata(cfg, alpha_source_timecourses_hilbert_dummy_bis);
   alpha_source_timecourses_hilbert_downsampled_bis.label = alpha_source_timecourses_hilbert.label; % we don't need dummy channels anymore
   
   %%% Not saved for disk space reasons
 %   save([subjectdata.subjectdir filesep 'HC_and_LC_' name_band_freq '_source_timecourses_hilbert_downsampled_bis'], 'alpha_source_timecourses_hilbert_downsampled_bis', '-v7.3')
 
 
      %% Rejection of the bad segments of the continuous signal
     disp('Rejection bad segments...');

 %Calculate the lengh of each run of continuous data
 for S = [1 2]  %Pour les 2 sessions
         datadir = ['subjectdata.datadir' int2str(S)];
         for R = [1 2]  %Pour les 2 runs
             if i == 1 || i == 6 || i == 10 || i == 11 || i == 13 || i == 15 || i == 16 || i == 20 || i == 22;
                  input_filename = [subjectdata.subjectdir filesep eval(datadir) filesep 'run0' int2str(R) '_tsss_realign_ica_eye_corrected_grad_ica_cardio_corrected_grad' '_continuous_' name_band_freq '.mat']
             else
                 input_filename = [subjectdata.subjectdir filesep eval(datadir) filesep 'run0' int2str(R) '_tsss_ica_eye_corrected_grad_ica_cardio_corrected_grad' '_continuous_' name_band_freq '.mat']
             end
             
             condition = eval(['subjectdata.condition' int2str(S) int2str(R)]);
             load(input_filename)
             
             if strcmp(condition,'High_contact') & R == 1
                 duration_HC_R1 = data_continuous_filtered.sampleinfo(1,2);
                 
             elseif strcmp(condition,'High_contact') & R == 2
                 duration_HC_R2 = data_continuous_filtered.sampleinfo(1,2);
                 
             elseif strcmp(condition,'Low_contact') & R == 1
                 duration_LC_R1 = data_continuous_filtered.sampleinfo(1,2);
                 
             elseif strcmp(condition,'Low_contact') & R == 2
                 duration_LC_R2 = data_continuous_filtered.sampleinfo(1,2);
                 
             end
             clear data_continuous_filtered 
         end
     end
     
 VERIF = duration_HC_R1 + duration_HC_R2 + duration_LC_R1 + duration_LC_R2;
 %------------------------------------------------------------------------------------------------------------------------------------------------------------------
  %Extract the good segments of data   
  %load([subjectdata.subjectdir filesep 'HC_and_LC_' name_band_freq '_source_timecourses_hilbert_downsampled_bis'])
  nsecs_artf = 2;
   for S = [1 2]  %Pour les 2 sessions
         	for R = [1 2]  %Pour les 2 runs   
                      
                     condition = eval(['subjectdata.condition' int2str(S) int2str(R)]);
                    % ---------------------------------------------------
                     if strcmp(condition,'High_contact') & R == 1
                     data_no_muscle_HC_R1 = eval(['subjectdata.data_no_muscle' int2str(S) int2str(R)]);
 
                     data_no_muscle_HC_R1(1,1) = data_no_muscle_HC_R1(1,1) + 1000;  %start 1second after the beginning of the session
                     
                     Nb_segment_HC_R1 = length (data_no_muscle_HC_R1(:, 1));
                     data_no_muscle_HC_R1(Nb_segment_HC_R1, 2) = (data_no_muscle_HC_R1(Nb_segment_HC_R1, 2)) - 1000; %cut 5seconds before the end of the session
                       
                     
                     for n = 1 : Nb_segment_HC_R1
                              
                             clean_HC_R1.trial{n} = alpha_source_timecourses_hilbert_downsampled_bis.trial{1, 1}(:, (round((data_no_muscle_HC_R1(n,1))/1000)+nsecs_artf):(round((data_no_muscle_HC_R1(n,2))/1000)-nsecs_artf));
                     
                             for x = 1 : length (alpha_source_timecourses_hilbert_downsampled_bis.mom)
                                 debut = 1;
                                 fin = length(alpha_source_timecourses_hilbert_downsampled_bis.mom{1, 1}(1,((round((data_no_muscle_HC_R1(n,1))/1000)+nsecs_artf):(round((data_no_muscle_HC_R1(n,2))/1000)-nsecs_artf))));
                                 clean_HC_R1.mom{x, n}(1,(debut:fin))   = alpha_source_timecourses_hilbert_downsampled_bis.mom{x, 1}(1,((round((data_no_muscle_HC_R1(n,1))/1000)+nsecs_artf):(round((data_no_muscle_HC_R1(n,2))/1000)-nsecs_artf)));
                                 clean_HC_R1.mom{x, n}(2,(debut:fin))   = alpha_source_timecourses_hilbert_downsampled_bis.mom{x, 1}(2,((round((data_no_muscle_HC_R1(n,1))/1000)+nsecs_artf):(round((data_no_muscle_HC_R1(n,2))/1000)-nsecs_artf)));
                                 clean_HC_R1.mom{x, n}(3,(debut:fin))   = alpha_source_timecourses_hilbert_downsampled_bis.mom{x, 1}(3,((round((data_no_muscle_HC_R1(n,1))/1000)+nsecs_artf):(round((data_no_muscle_HC_R1(n,2))/1000)-nsecs_artf)));
                             end
                     end
                         
                      
                       conca_clean_HC_R1.trial = {[clean_HC_R1.trial{:}]};  
                       
                       for  x = 1 : length (clean_HC_R1.mom)
                           conca_clean_HC_R1.mom{x,1} = [clean_HC_R1.mom{x, :}];
                       end
                       
                      clear data_no_muscle_HC_R1 Nb_segment_HC_R1 start_HC_R1 clean_HC_R1 
                    % ---------------------------------------------------
                     elseif strcmp(condition,'High_contact') & R == 2
                     data_no_muscle_HC_R2 = eval(['subjectdata.data_no_muscle' int2str(S) int2str(R)]);
 
                     data_no_muscle_HC_R2 = data_no_muscle_HC_R2 + duration_HC_R1;
                     data_no_muscle_HC_R2(1,1) = data_no_muscle_HC_R2(1,1) + 1000 ;
                     
                     Nb_segment_HC_R2 = length (data_no_muscle_HC_R2(:, 1));
                     data_no_muscle_HC_R2(Nb_segment_HC_R2, 2) = (data_no_muscle_HC_R2(Nb_segment_HC_R2, 2)) - 1000;
                     
                     for n = 1 : Nb_segment_HC_R2
                              
                              clean_HC_R2.trial{n} = alpha_source_timecourses_hilbert_downsampled_bis.trial{1, 1}(:, (round((data_no_muscle_HC_R2(n,1))/1000)+nsecs_artf):(round((data_no_muscle_HC_R2(n,2))/1000)-nsecs_artf)); 
                     
                             for x = 1 : length (alpha_source_timecourses_hilbert_downsampled_bis.mom)
                             debut = 1;
                             fin = length(alpha_source_timecourses_hilbert_downsampled_bis.mom{1, 1}(1,((round((data_no_muscle_HC_R2(n,1))/1000)+nsecs_artf):(round((data_no_muscle_HC_R2(n,2))/1000)-nsecs_artf))));
                             clean_HC_R2.mom{x, n}(1,(debut:fin))   = alpha_source_timecourses_hilbert_downsampled_bis.mom{x, 1}(1,((round((data_no_muscle_HC_R2(n,1))/1000)+nsecs_artf):(round((data_no_muscle_HC_R2(n,2))/1000)-nsecs_artf)));
                             clean_HC_R2.mom{x, n}(2,(debut:fin))   = alpha_source_timecourses_hilbert_downsampled_bis.mom{x, 1}(2,((round((data_no_muscle_HC_R2(n,1))/1000)+nsecs_artf):(round((data_no_muscle_HC_R2(n,2))/1000)-nsecs_artf)));
                             clean_HC_R2.mom{x, n}(3,(debut:fin))   = alpha_source_timecourses_hilbert_downsampled_bis.mom{x, 1}(3,((round((data_no_muscle_HC_R2(n,1))/1000)+nsecs_artf):(round((data_no_muscle_HC_R2(n,2))/1000)-nsecs_artf)));
                             end      
                     end
                         
                       
                       conca_clean_HC_R2.trial = {[clean_HC_R2.trial{:}]} ; 
                       
                       for  x = 1 : length (clean_HC_R2.mom)
                       conca_clean_HC_R2.mom{x,1} = [clean_HC_R2.mom{x, :}];
                       end 
                         
                      clear data_no_muscle_HC_R2 Nb_segment_HC_R2 start_HC_R2 clean_HC_R2
                   %  ---------------------------------------------------
                     elseif strcmp(condition,'Low_contact') & R == 1
                      data_no_muscle_LC_R1 = eval(['subjectdata.data_no_muscle' int2str(S) int2str(R)]);
                      
                     data_no_muscle_LC_R1 = data_no_muscle_LC_R1 + duration_HC_R2 + duration_HC_R1;
                     
                     data_no_muscle_LC_R1(1,1) = data_no_muscle_LC_R1(1,1) + 1000;  
                     
                     Nb_segment_LC_R1 = length (data_no_muscle_LC_R1(:, 1));
                     data_no_muscle_LC_R1(Nb_segment_LC_R1, 2) = (data_no_muscle_LC_R1(Nb_segment_LC_R1, 2)) - 1000;
                     
                     for n = 1 : Nb_segment_LC_R1
                              
                             clean_LC_R1.trial{n} = alpha_source_timecourses_hilbert_downsampled_bis.trial{1, 1}(:, (round((data_no_muscle_LC_R1(n,1))/1000)+nsecs_artf):(round((data_no_muscle_LC_R1(n,2))/1000)-nsecs_artf)); 
                     
                             for x = 1 : length (alpha_source_timecourses_hilbert_downsampled_bis.mom)
                             debut = 1;
                             fin = length(alpha_source_timecourses_hilbert_downsampled_bis.mom{1, 1}(1,((round((data_no_muscle_LC_R1(n,1))/1000)+nsecs_artf):(round((data_no_muscle_LC_R1(n,2))/1000)-nsecs_artf))));
                             clean_LC_R1.mom{x, n}(1,(debut:fin))   = alpha_source_timecourses_hilbert_downsampled_bis.mom{x, 1}(1,((round((data_no_muscle_LC_R1(n,1))/1000)+nsecs_artf):(round((data_no_muscle_LC_R1(n,2))/1000)-nsecs_artf)));
                             clean_LC_R1.mom{x, n}(2,(debut:fin))   = alpha_source_timecourses_hilbert_downsampled_bis.mom{x, 1}(2,((round((data_no_muscle_LC_R1(n,1))/1000)+nsecs_artf):(round((data_no_muscle_LC_R1(n,2))/1000)-nsecs_artf)));
                             clean_LC_R1.mom{x, n}(3,(debut:fin))   = alpha_source_timecourses_hilbert_downsampled_bis.mom{x, 1}(3,((round((data_no_muscle_LC_R1(n,1))/1000)+nsecs_artf):(round((data_no_muscle_LC_R1(n,2))/1000)-nsecs_artf)));
                             end      
                     end
                         
                       
                       conca_clean_LC_R1.trial = {[clean_LC_R1.trial{:}]};  
                       
                       for  x = 1 : length (clean_LC_R1.mom)
                       conca_clean_LC_R1.mom{x,1} = [clean_LC_R1.mom{x, :}];
                       end 
                       
                      clear data_no_muscle_LC_R1 Nb_segment_LC_R1 start_LC_R1 clean_LC_R1
                    % ---------------------------------------------------
                     elseif strcmp(condition,'Low_contact') & R == 2
                     data_no_muscle_LC_R2 = eval(['subjectdata.data_no_muscle' int2str(S) int2str(R)]);
                     
                     data_no_muscle_LC_R2 = data_no_muscle_LC_R2 + duration_HC_R2 + duration_HC_R1 + duration_LC_R1;
                     data_no_muscle_LC_R2(1,1) = data_no_muscle_LC_R2(1,1) + 1000 ;
                     
                     Nb_segment_LC_R2 = length (data_no_muscle_LC_R2(:, 1));
                     data_no_muscle_LC_R2(Nb_segment_LC_R2, 2) = (data_no_muscle_LC_R2(Nb_segment_LC_R2, 2)) - 1000; 
                     
                     for n = 1 : Nb_segment_LC_R2
                              
                              clean_LC_R2.trial{n} = alpha_source_timecourses_hilbert_downsampled_bis.trial{1, 1}(:, (round((data_no_muscle_LC_R2(n,1))/1000)+nsecs_artf):(round((data_no_muscle_LC_R2(n,2))/1000)-nsecs_artf)); 
                     
                             for x = 1 : length (alpha_source_timecourses_hilbert_downsampled_bis.mom)
                             debut = 1;
                             fin = length(alpha_source_timecourses_hilbert_downsampled_bis.mom{1, 1}(1,((round((data_no_muscle_LC_R2(n,1))/1000)+nsecs_artf):(round((data_no_muscle_LC_R2(n,2))/1000)-nsecs_artf))));
                             clean_LC_R2.mom{x, n}(1,(debut:fin))   = alpha_source_timecourses_hilbert_downsampled_bis.mom{x, 1}(1,((round((data_no_muscle_LC_R2(n,1))/1000)+nsecs_artf):(round((data_no_muscle_LC_R2(n,2))/1000)-nsecs_artf)));
                             clean_LC_R2.mom{x, n}(2,(debut:fin))   = alpha_source_timecourses_hilbert_downsampled_bis.mom{x, 1}(2,((round((data_no_muscle_LC_R2(n,1))/1000)+nsecs_artf):(round((data_no_muscle_LC_R2(n,2))/1000)-nsecs_artf)));
                             clean_LC_R2.mom{x, n}(3,(debut:fin))   = alpha_source_timecourses_hilbert_downsampled_bis.mom{x, 1}(3,((round((data_no_muscle_LC_R2(n,1))/1000)+nsecs_artf):(round((data_no_muscle_LC_R2(n,2))/1000)-nsecs_artf)));
                             end      
                     end
                         
                       conca_clean_LC_R2.trial = {[clean_LC_R2.trial{:}]};  
                       
                       for  x = 1 : length (clean_LC_R2.mom)
                       conca_clean_LC_R2.mom{x,1} = [clean_LC_R2.mom{x, :}];
                       end 
                         
                      clear data_no_muscle_LC_R2 Nb_segment_LC_R2 start_LC_R2 clean_LC_R2
                    % ---------------------------------------------------
                     
                     end          
             end
   end
   VERIF_bis = size(alpha_source_timecourses_hilbert_downsampled_bis.trial{1, 1});
   
   % Concatenate data
   concat_data = [conca_clean_HC_R1.trial{1} conca_clean_HC_R2.trial{1} conca_clean_LC_R1.trial{1} conca_clean_LC_R2.trial{1}];
   
   % Check if outliers remaining (outliers identified in the average of all voxels, >avg+4std)
   nnz_data = find(concat_data(:, 1) ~= 0);
   avg_concat_data = mean(concat_data(nnz_data, :), 1);
   gd_avg_concat_data = mean(avg_concat_data);
   gd_std_concat_data = std(avg_concat_data);
   outl_pts = find(avg_concat_data > gd_avg_concat_data + 4*gd_std_concat_data);
   concat_data(:, outl_pts) = [];
   
   % Create a fieldtrip-like data structure
   alpha_source_timecourses_hilbert_clean_bis.trial{1}  = concat_data;
   alpha_source_timecourses_hilbert_clean_bis.time   = {[1:1:(length(alpha_source_timecourses_hilbert_clean_bis.trial{1,1}(1,:)))]}; 
   alpha_source_timecourses_hilbert_clean_bis.pos    = alpha_source_timecourses_hilbert_downsampled_bis.pos;
   alpha_source_timecourses_hilbert_clean_bis.filter = alpha_source_timecourses_hilbert_downsampled_bis.filter;
   alpha_source_timecourses_hilbert_clean_bis.label  = alpha_source_timecourses_hilbert_downsampled_bis.label;
   
   for v = 1 : length (alpha_source_timecourses_hilbert_downsampled_bis.mom)
       alpha_source_timecourses_hilbert_clean_bis.mom {v, 1}    = [conca_clean_HC_R1.mom{v, 1} conca_clean_HC_R2.mom{v, 1} conca_clean_LC_R1.mom{v, 1} conca_clean_LC_R2.mom{v, 1}];
       alpha_source_timecourses_hilbert_clean_bis.mom {v, 1}(:, outl_pts) = [];
   end
   
   
 % Add conditions informations to the concatenated data
 n1 = length(conca_clean_HC_R1.trial{1}(1,:));
 n2 = length(conca_clean_HC_R2.trial{1}(1,:));
 n3 = length(conca_clean_LC_R1.trial{1}(1,:));
 n4 = length(conca_clean_LC_R2.trial{1}(1,:));
 
 alpha_source_timecourses_hilbert_clean_bis.condition = [1: (n1 + n2 + n3+ n4)];
 alpha_source_timecourses_hilbert_clean_bis.condition(1, (1:n1)) = 11;
 alpha_source_timecourses_hilbert_clean_bis.condition(1, ((1+n1):(n1+n2))) = 12;
 alpha_source_timecourses_hilbert_clean_bis.condition(1, ((1+n1+n2):(n1+n2+n3))) = 21;
 alpha_source_timecourses_hilbert_clean_bis.condition(1, ((1+n1+n2+n3):(n1+n2+n3+n4))) = 22;
 alpha_source_timecourses_hilbert_clean_bis.condition(outl_pts) = [];
 
 %%% Not saved for disk space reasons
 % save([subjectdata.subjectdir filesep 'HC_and_LC_' name_band_freq '_source_timecourses_hilbert_clean_bis'], 'alpha_source_timecourses_hilbert_clean_bis', '-v7.3')
   
 % Plot raw and clean data
 disp('Plotting raw data...');
 figure('color', [1 1 1]);
 subplot(2, 1, 1);
 nnz_data = find(alpha_source_timecourses_hilbert_downsampled_bis.trial{1, 1}(:, 1) ~= 0);
 avg_raw = mean(alpha_source_timecourses_hilbert_downsampled_bis.trial{1, 1}(nnz_data, :));
 plot(avg_raw);
 title(['Subj ' num2str(i) ', n pts=' num2str(length(avg_raw))]);
 xlim([0 length(avg_raw)]);
 subplot(2, 1, 2);
 avg_concat_data = mean(concat_data(nnz_data, :), 1);
 plot(avg_concat_data);
 title(['Subj ' num2str(i) ', n pts=' num2str(length(avg_concat_data)) ', n outl=' num2str(length(outl_pts)) ' - clean']);
 xlim([0 length(avg_raw)]);
 print('-dpng', ['~/../datalinks/SELFI/Networks_Group_comparison/' name_band_freq '/visu/Timecourse_AvgVox_RawClean/Subj' num2str(i)]);

 
 % Clear variables
 clear conca_clean_HC_R1 conca_clean_HC_R2 conca_clean_LC_R1 conca_clean_LC_R2
 
 
end
