%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Batch_MRI_preprocessing_and_head_model
%%%
%%%  For each subject preprocess the MRI: units conversion, reslicing, realignement and segmentation
%%%  and create a head model
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Individual MRI Preprocessing (realignement in MEG coordinates system) : interactive method

clear all
close all

subjects = [1, 2, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 19, 20, 21, 22]; 
        
for i = subjects
  try
        eval(['subject' int2str(i)]);
        Outputdir = eval('subjectdata.subjectdir')
        mkdir(eval('Outputdir'),'Anatomy')      
        Raw_MRI_file = [ '/pclxserver/raid11/data/SELFI/RAW_MRI' filesep eval('subjectdata.IRMdate') '_SELFI_S' eval('subjectdata.subjectnr') filesep 'S02_t1_mpr_sag_0_8iso' filesep 's_S02_t1_mpr_sag_0_8iso.nii.gz']
        
        % Load the MRI
        mri = ft_read_mri(Raw_MRI_file, 'dataformat', 'nifti'); %Nifti is in mm / Neuromag in m / CTF in cm / see http://www.fieldtriptoolbox.org/faq/how_are_the_different_head_and_mri_coordinate_systems_defined
       
        % Realign the MRI into Neuromag by help of fiducial points
        % Making sure you know which side is the right side (ours MRI are
        % inversed, the left is on the right, radiological convention)
        % assign the nasion (pressing "n"), left ("l") and right ("r") with the crosshairs on
        % the ear markers. Then finish with "q".
        cfg = [];
        cfg.method = 'interactive';
        cfg.coordsys = 'neuromag';
        mri_realigned = ft_volumerealign(cfg,mri);  
        
        save([eval('Outputdir') filesep 'Anatomy/mri_realigned.mat'], 'mri_realigned')
        
  catch
      disp(['Something was wrong with Subject' int2str(i) 'for Realignement! Continuing with next in line']);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Individual MRI Preprocessing (units conversion, reslicing, realignement and segmentation)

clear all
close all 

subjects = [1, 2, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 19, 20, 21, 22]; 


for i = subjects
  try
        eval(['subject' int2str(i)]);
        Outputdir = eval('subjectdata.subjectdir')
       
        % Load the realigned MRI
        mri_realigned_filename = [eval('Outputdir') filesep 'Anatomy/mri_realigned.mat']
        load(mri_realigned_filename)

        % Reslice the MRI  
        % i.e. to interpolate the anatomy onto a new 3D grid that is aligned with the axes of the coordinate system. 
        % This prevents problems such as describe here :
        % http://www.fieldtriptoolbox.org/faq/why_does_my_anatomical_mri_show_upside-down_when_plotting_it_with_ft_sourceplot
        % and http://www.fieldtriptoolbox.org/faq/my_mri_is_upside_down_is_this_a_problem
        cfg = [];
        mri_realigned_resliced = ft_volumereslice(cfg, mri_realigned)
        
        save([ eval('Outputdir') filesep 'Anatomy/mri_realigned_resliced.mat'], 'mri_realigned_resliced')
        
        % converting from 'mm' to 'm' units (car neuromag en m)
        mri_realigned_resliced_meters = ft_convert_units(mri_realigned_resliced, 'm');
        
        save([ eval('Outputdir') filesep 'Anatomy/mri_realigned_resliced_meters.mat'], 'mri_realigned_resliced_meters')

      %% Segmentation of the preprocessed MRI
      
      % Segment the realigned and resliced MRI (en m) - took 153 seconds
      tic
      cfg = [];
      cfg.output    = {'brain'};
      [segmentedmri] = ft_volumesegment(cfg, mri_realigned_resliced_meters); 
      
      save( [eval('Outputdir') filesep 'Anatomy/segmentedmri.mat'], 'segmentedmri')
      
      t = toc;
      fprintf(['\n###############\nFinished - segmentation took ' num2str(t) ' seconds.\n\n'])
      
      
      
%       % ! Make a figure of the mri and segmented volumes
%       segmentedmri           = segmentedmri;
%       segmentedmri.transform = mri_realigned_resliced_meters.transform;
%       segmentedmri.anatomy   = mri_realigned_resliced_meters.anatomy;  
%       cfg                    = [];
%       cfg.funparameter       = 'brain';
%       ft_sourceplot(cfg, segmentedmri);

      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
      %% Create the head model from the individual MRI (needed by the beamforming to compute a forward model)

      % Create the head model from the individual MRI
      cfg         = [];
      cfg.method  = 'singleshell';
      headmodel   = ft_prepare_headmodel(cfg, segmentedmri);
      headmodel   = ft_convert_units( headmodel, 'm'); % mm to m, since the grid will also be expressed in m
      
      save([eval('Outputdir') filesep 'Anatomy/headmodel.mat'],'headmodel')
      


