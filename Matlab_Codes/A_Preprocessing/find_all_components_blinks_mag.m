function [components_to_remove, coeff] = find_all_components_blinks_mag(input_filename, data_no_muscle)


cfg         = [];
cfg.dataset = input_filename ;
cfg.channel = {'BIO002'} ;
% cfg.hpfilter = 'yes';
% cfg.hpfreq = 0.1;
% cfg.hpfiltord = 3;
% cfg.lpfilter = 'yes';
% cfg.lpfreq = 30;
% cfg.lpfiltord = 3;

data        = ft_preprocessing(cfg);

% data.trial{1} = [data.trial{1}(:,1:10000) data.trial{1}(:,11000:15000) data.trial{1}(:,200000:end)];

Nb_segment = length (data_no_muscle(:,1))
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

%% downsample BIO channels 
cfg                  = [];
cfg.resamplefs       = 100;
cfg.detrend          = 'no';
data_downsample      = ft_resampledata(cfg, data);
data = data_downsample;


%% load components file
[path,name,~] = fileparts(input_filename) ;
components_filename = [ path '/' name '_ica_comp_blinks.mat' ];
load (components_filename);


%% compute the correlation coef between BIO channels and each components

components_to_remove = [];
k=0;
coeff = NaN(length(comp.trial{1}(:,1)),2); %create matrice de Not-a-Number, avec nbr ligne = nbr ICA composante, nbr colonne = 2 (une colonne avec coeff correlation, une colonne binaire avec 1 si corr> 0.5)

for i = 1 : length(comp.trial{1}(:,1)) %nbr de composantes
    
    C1 = corrcoef(comp.trial{1}(i,:), data.trial{1}(1,:));     % corr entre composante nbr i (the all time points) x channel Bio (the all time points), prend la forme d'une matrice de corr 2*2 avec (avec valeur 1 dans la diagonale) 
    coeff(i,1) = C1(1,2);    %recup dans la matrice de NaN "coeff" (de taille 69*2) sur ligne i dans la 1?re colonne, le coeff de corr
    if abs(coeff(i,1)) > 0.5
        coeff(i,2) = 1;
        k=k+1;
        components_to_remove(k)=i;
    end
    
 %% write the output

[path,name,] = fileparts(input_filename) ;
components_to_remove_filename = [ path '/' name '_components_to_remove_blinks.mat' ] ;
display( [ 'Saving ' components_to_remove_filename ] )
save( components_to_remove_filename, 'components_to_remove' );   
end

components_to_remove