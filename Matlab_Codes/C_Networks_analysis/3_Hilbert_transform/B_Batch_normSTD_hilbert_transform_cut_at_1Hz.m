%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Batch_normSTD_hilbert_transform_cut_at_1Hz
%%%
%%% In detail, it takes the source time courses from subject i 
%%% concatenated from the 4 runs (computed in Batch_concatenation4runs)  
%%%   1) load the downsampled Hilbert envelopes
%%%   2) normalizes them
%%%   3) and smooth them
%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

subjects = [1, 2, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 19, 20, 21, 22]  %pclx41

name_band_freq = 'alpha'

for i = subjects
    
 eval(['subject' int2str(i)])
 Outputdir = eval('subjectdata.subjectdir')
  
load([subjectdata.subjectdir filesep 'HC_and_LC_' name_band_freq '_source_timecourses_hilbert_clean_bis'])

%---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  % Normalization
  alpha_source_timecourses_hilbert_downsampled_norm_bis = alpha_source_timecourses_hilbert_clean_bis;
    for ngrid = 1:size(alpha_source_timecourses_hilbert_downsampled_norm_bis.trial{1},1)
        alpha_source_timecourses_hilbert_downsampled_norm_bis.trial{1}(ngrid,:) = (alpha_source_timecourses_hilbert_downsampled_norm_bis.trial{1}(ngrid,:) - mean(alpha_source_timecourses_hilbert_downsampled_norm_bis.trial{1}(ngrid,:))) ...
                                                                     / std(alpha_source_timecourses_hilbert_downsampled_norm_bis.trial{1}(ngrid,:));                                                           
    end
        
    save([subjectdata.subjectdir filesep 'HC_and_LC_' name_band_freq '_source_timecourses_hilbert_downsampled_norm_bis'], 'alpha_source_timecourses_hilbert_downsampled_norm_bis', '-v7.3')

  
  clear alpha_source_timecourses_hilbert_clean_bis
  
 %----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  % Load the grid and the hilbert-transformed-source-timecourse
  load ([Outputdir filesep 'Anatomy/grid1cm.mat'],'grid1cm')
   
    
        alpha_source_timecourses_hilbert_downsampled_norm_smoothed = alpha_source_timecourses_hilbert_downsampled_norm_bis;

        % smooth data timestep by timestep
        for timestep = 1:size(alpha_source_timecourses_hilbert_downsampled_norm_smoothed.trial{1},2)
            
            data = alpha_source_timecourses_hilbert_downsampled_norm_smoothed.trial{1}(:,timestep);

            % reshape data into a 3D matrix with the dimensions of the grid
            data_reshape =[];

            % wrap every 3-tupel in cell first
            for ngrid=1:length(grid1cm.pos)
                data_reshape = [data_reshape; {data(ngrid,:)}];
            end
            
            data_reshape = reshape(data,[grid1cm.dim(1), grid1cm.dim(2), grid1cm.dim(3)]);

            % set NaN values to zero to avoid problems with imgaussfilt3
            data_reshape(isnan(data_reshape)) =0;

            % smooth the data
            data_smooth = imgaussfilt3(data_reshape, [1 1 1],'FilterDomain', 'spatial');

            % reshape back into original dimensions
            data_smooth = reshape(data_smooth,[length(grid1cm.pos) 1 1]);
            alpha_source_timecourses_hilbert_downsampled_norm_smoothed.trial{1}(:,timestep) = data_smooth;

        end
        
        save([subjectdata.subjectdir filesep 'HC_and_LC_' name_band_freq '_source_timecourses_hilbert_downsampled_norm_smoothed'], 'alpha_source_timecourses_hilbert_downsampled_norm_smoothed', '-v7.3')
        clear grid1cm alpha_source_timecourses_hilbert_downsampled_norm_bis alpha_source_timecourses_hilbert_downsampled_norm_smoothed data_reshape ngrid data_smooth data timestep
 
  
end



