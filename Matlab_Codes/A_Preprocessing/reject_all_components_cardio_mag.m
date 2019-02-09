function [ data_corrected, output_filename ] = reject_all_components_cardio( input_filename, components_filename, components_to_remove )

%[~,name,~] = fileparts(input_filename) ;

[path,name,~] = fileparts(input_filename) ;
components_filename = [ path '/' name '_ica_comp_cardio.mat' ];

%% preprocessing of example dataset
cfg         = [];
cfg.dataset = input_filename ;
%cfg.channel = 'MEG*1';
data        = ft_preprocessing(cfg);

%% read the ICA components


components            = load( components_filename ) ;
load( components_to_remove ) 

% correction
cfg                           = [];
cfg.component                 = components_to_remove ; % to be removed component(s)
data_corrected                = ft_rejectcomponent(cfg, components.comp, data );



% reordering the channels
MEG0111_id=-1;
for i=1:length(data.label)
    if strcmp( data.label{i}, 'MEG0111' )
        MEG0111_id=i ;
        break ;
    end
end
if MEG0111_id == -1
    error( 'MEG0111 not found' )
end
MEG2643_id=-1 ;
for i=MEG0111_id+1:length(data.label)
    if strcmp( data.label{i}, 'MEG2643' )
        MEG2643_id=i ;
        break ;
    end
end
if MEG2643_id == -1
    error( 'MEG2643 not found' )
end


cfg         = [];
cfg.channel = data.label(1:MEG0111_id - 1);
data_corrected_1 = ft_preprocessing(cfg, data_corrected);

cfg.channel = data.label(MEG0111_id:MEG2643_id);
data_corrected_2 = ft_preprocessing(cfg, data_corrected);

cfg.channel = data.label(MEG2643_id+1:length(data.label));
data_corrected_3 = ft_preprocessing(cfg, data_corrected);

cfg         = [];
data_corrected = ft_appenddata(cfg, data_corrected_1, data_corrected_2, data_corrected_3) ;
data_corrected.hdr = data.hdr ;

% writing the results
[path,name,~] = fileparts(input_filename) ;
 
output_filename = [ path '/' name '_ica_cardio_corrected' ];


fieldtrip2fiff( [ output_filename '.fif' ], data_corrected );
