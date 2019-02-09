function [ n_comp, comp, components_filename ] = do_ica_cardio_mag( input_filename, data_no_muscle )


%% load the data
cfg         = [];
cfg.dataset = input_filename ;
cfg.channel = 'MEG*1';
%cfg.hpfilter = 'yes';
%cfg.hpfreq = 0.1;
%cfg.hpfiltord = 3;
% cfg.lpfilter = 'yes';
% cfg.lpfreq = 30;
% cfg.lpfiltord = 3;

data        = ft_preprocessing(cfg);


Nb_segment = length (data_no_muscle(:,1)) % nbr de segments de data propre
Nb_sample_tot = 0

for i = 1 : Nb_segment
Nb_sample_i = (data_no_muscle(i,2) - data_no_muscle(i,1)) +1
i = i+1
Nb_sample_tot = Nb_sample_tot + Nb_sample_i 
end

cleandata.trial{1} = NaN(length(data.trial{1}(:,1)),Nb_sample_tot);
end_column = 0

for i = 1 : Nb_segment
    start_sample = data_no_muscle(i,1) %start of the sample to add
    end_sample = data_no_muscle(i,2) %end of the sample to add
   
    segment = [data.trial{1}(:,((start_sample):(end_sample)))] %samples to add
    
    k= (end_sample - start_sample) %nombre de samples to add
    start_column = end_column + 1
    end_column = start_column + k 
    cleandata.trial{1}(:, ((start_column):(end_column))) = [segment]
    i = i+1   
end   

data.trial{1} = cleandata.trial{1}
data.sampleinfo (1, 2) = Nb_sample_tot

cleandata.time{1} = NaN(1,Nb_sample_tot);

for i = 1 : Nb_sample_tot
    cleandata.time{1}(1, i) = data.time{1}(1, i)
    i = i+1   
end   

data.time{1} = cleandata.time{1}

%% downsample the data to speed up the next step
cfg                  = [];
cfg.resamplefs       = 100;
cfg.detrend          = 'no';
data_downsample      = ft_resampledata(cfg, data);

%% add sampleinfo to use ft_databrowser
data_downsample      = ft_checkdata( data_downsample, 'hassampleinfo', 'yes') ;

%% run the ICA
cfg              = [];
cfg.method       = 'runica';
n_comp           = rank(squeeze(data_downsample.trial{1}) * squeeze(data_downsample.trial{1})');
cfg.numcomponent = n_comp;
comp            = ft_componentanalysis(cfg, data_downsample);


%% write the output
[path,name,] = fileparts(input_filename) ;
components_filename = [ path '/' name '_ica_comp_cardio.mat' ] ;
display( [ 'Saving ' components_filename ] )
save( components_filename, 'comp' );

end