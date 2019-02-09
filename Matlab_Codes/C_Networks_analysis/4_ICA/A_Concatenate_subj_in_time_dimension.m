%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SCRIPT_Concatenate_subj_in_time_dimension
%%%
%%% iteratively load the hilbert source time courses for each subject
%%% containing the 4 runs, and appends them
%%% into a whole vector across subjects across conditions. The order of the
%%% concatenated data is then the following:
%%% [Subject1_HCR1 Subject1_HCR2 Subject1_LCR1 Subject1_LCR2
%%% Subject2_HCR1....]
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
name_band_freq = 'theta';

id_subjects = [1, 2, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 19, 20, 21, 22];

Output = ['/pclxserver/raid11/data/SELFI/Networks_Group_comparison/' name_band_freq '/SourceData'];

load(['/pclxserver/raid11/data/SELFI/selfi01_s01' filesep 'HC_and_LC_' name_band_freq '_source_timecourses_hilbert_downsampled_norm_smoothed'])

alpha_hilbert_across_subjects_across_groups = [];
alpha_hilbert_across_subjects_across_groups.label = alpha_source_timecourses_hilbert_downsampled_norm_smoothed.label;
alpha_hilbert_across_subjects_across_groups.pos = alpha_source_timecourses_hilbert_downsampled_norm_smoothed.pos;
alpha_hilbert_across_subjects_across_groups.trial = alpha_source_timecourses_hilbert_downsampled_norm_smoothed.trial;
alpha_hilbert_across_subjects_across_groups.condition = alpha_source_timecourses_hilbert_downsampled_norm_smoothed.condition;

for i = 2:length(id_subjects)
    
    subjects = id_subjects(i);
    eval(['subject' int2str(subjects)])
    
    Inputdir = eval('subjectdata.subjectdir');
    load([subjectdata.subjectdir filesep 'HC_and_LC_' name_band_freq '_source_timecourses_hilbert_downsampled_norm_smoothed'])
    
    alpha_hilbert_across_subjects_across_groups.trial = {[alpha_hilbert_across_subjects_across_groups.trial{1} alpha_source_timecourses_hilbert_downsampled_norm_smoothed.trial{1}]};
    alpha_hilbert_across_subjects_across_groups.condition = [alpha_hilbert_across_subjects_across_groups.condition alpha_source_timecourses_hilbert_downsampled_norm_smoothed.condition];
 
end

alpha_hilbert_across_subjects_across_groups.time = {[1:1:size(alpha_hilbert_across_subjects_across_groups.trial{1},2)]};

% fill empty entries with zeros and fill them with zeros (to avoid problems
% when running PCA/ICA)
[rowsnan, ~] = find(isnan(alpha_hilbert_across_subjects_across_groups.trial{1}));
alpha_hilbert_across_subjects_across_groups.trial{1}(rowsnan,:) = 0;

% Rename variable to have a neutral name regarding the frequency band
data_hilbert_all = alpha_hilbert_across_subjects_across_groups;

% Save
save([Output filesep name_band_freq '_hilbert_across_subjects_across_groups_smooth'], 'data_hilbert_all', '-v7.3')


