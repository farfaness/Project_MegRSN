%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SCRIPT_Visu
%%% 
%%% Script to visualise the components of the ICA to correct blinks and cardio for magnometers
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Start by defining the inputs:

% subject_number= 2 ; session_number= 2; run_number = 2

function [] = visualize_ica_components(subject_number, session_number, run_number)

%%load components file

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


input_filename = [subjectdata.subjectdir filesep datadir filesep run '_tsss.fif']
[path,name,~] = fileparts(input_filename) ;
components_filename = [ path '/' name '_ica_comp.mat' ];
load (components_filename);

%% plot the components for visual inspection figure
cfg = [];
cfg.component = [1:10];       % specify the component(s) that should be plotted
cfg.colormap =     jet           % See COLORMAP
cfg.layout = '/lena13/home_users/users/hazem/fieldtrip-20161231/template/layout/neuromag306mag.lay';
cfg.comment   = 'no';
cfg.renderer = 'painters'
ft_topoplotIC(cfg, comp)

cfg = [];
cfg.layout = '/lena13/home_users/users/hazem/fieldtrip-20161231/template/layout/neuromag306mag.lay';
cfg.viewmode = 'component'
cfg.renderer = 'painters' 
ft_databrowser(cfg, comp)

%%Display the number of the ica components identyfied as artefacts
components_to_remove_filename = [ path '/' name '_components_to_remove.mat' ];
load (components_to_remove_filename); 
disp(['The artefact components are :' num2str(components_to_remove)])

end
