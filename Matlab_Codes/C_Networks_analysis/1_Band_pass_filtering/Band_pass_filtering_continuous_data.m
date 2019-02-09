function [] = Band_pass_filtering_continuous_data(input_filename, condition, lpfreq, hpfreq, name_band_freq )

%%%%% BAND_PASS_FILTERING 

%   Band_pass filter data in fonction of the specified band of frequency
%   computed on continuous data


%Determine condition : 1 for High Contact, 2 for Low Contact
if strcmp(condition,'High_contact')
     cond = 1;
elseif strcmp(condition,'Low_contact')
     cond = 2;
end

%Load the continous data
cfg         = [];
cfg.dataset = input_filename ;
data_continuous = ft_preprocessing(cfg);

%% compute band-pass filtering of the data
cfg = [];
cfg.lpfilter      = 'yes';
cfg.hpfilter      = 'yes';
cfg.lpfreq        = lpfreq;
cfg.hpfreq        = hpfreq;
[data_continuous_filtered] = ft_preprocessing(cfg, data_continuous);

%% Save the output 
[path,name,] = fileparts(input_filename) ;
Band_pass_filtered_filename = [ path '/' name '_continuous_' name_band_freq '.mat' ] ;
display( [ 'Saving ' Band_pass_filtered_filename  ] )
save( Band_pass_filtered_filename, 'data_continuous_filtered' );


end

