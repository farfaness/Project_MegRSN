%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%% Visualisation of power spectrum
%%%%% 
%%%%% 1)  Multiple channels power spectrum
%%%%% 2)  Topography on frequency bands : delta, theta, alpha, beta, gamma
%%%%% 3)  Power spectrum log10 on all the magnetometers averaged togethers
%%%%% 4)  Individual Power Spectrum for all the sujects and each runs (not averaged over channels)
%%%%% 5)  Individual Power Spectrum for all the subjects: COND*RUN (averaged over channels)
%%%%% 6)  Individual Power Spectrum for all the subjects: COND (averaged over runs and channels)
%%%%% 7)  Vizu Power Spectrum on the subject with order HL and conditions taken separatly (averaged over channels and runs)
%%%%% 8)  Vizu Power Spectrum on the subject with order LH and conditions taken separatly (averaged over channels and runs)
%%%%% 9)  Vizu Power Spectrum on the subject with order HL and conditions/runs taken separatly (averaged over channels)
%%%%% 10) Vizu Power Spectrum on the subject with order LH and conditions/runs taken separatly (averaged over channels)
%%%%% 11) Power Spectrum Cond*Order (average over subjects)
%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%% Multiple channels power spectrum
cfg = []
% cfg.zlim = [log10(-3e-27) log10(3e-27)]; %limite valeur Z (power)
cfg.ZScale = 'log'
cfg.showlabels   = 'yes';
cfg.layout = '/lena13/home_users/users/hazem/fieldtrip-20161231/template/layout/neuromag306mag.lay'   % for magnetometers only 
%cfg.parameter     = field to be plotted on y-axis (default depends on data.dimord) 'avg', 'powspctrm' or 'cohspctrm'
figure; ft_multiplotER(cfg, power_low_contact, power_high_contact);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Topography delta, theta, alpha, beta, gamma

% Topography delta
cfg = [];
cfg.xlim = [1 4] %[1 : 5 : 40];
cfg.zlim         = [-8e-28 8e-28];  %'maxmin' color for power
cfg.colormap =     jet ;
cfg.parameter    = 'powspctrm' 
%cfg.marker       = 'on';
cfg.comment = 'xlim';
cfg.commentpos = 'title';
cfg.layout = '/lena13/home_users/users/hazem/fieldtrip-20161231/template/layout/neuromag306mag.lay'
figure; ft_topoplotTFR(cfg, power_low_contact);

% Topography theta 
cfg = [];
cfg.xlim = [4 8] %[1 : 5 : 40];
cfg.zlim         = [-6e-28 6e-28];  %'maxmin' color for power
cfg.colormap =     jet ;
cfg.parameter    = 'powspctrm' 
%cfg.marker       = 'on';
cfg.comment = 'xlim';
cfg.commentpos = 'title';
cfg.layout = '/lena13/home_users/users/hazem/fieldtrip-20161231/template/layout/neuromag306mag.lay'
figure; ft_topoplotTFR(cfg, power_low_contact);

% Topography alpha 
cfg = [];
cfg.xlim = [8 13] %[1 : 5 : 40];
cfg.zlim         = [-9e-28 9e-28];  %'maxmin' color for power
cfg.colormap =     jet ;
cfg.parameter    = 'powspctrm' 
%cfg.marker       = 'on';
cfg.comment = 'xlim';
cfg.commentpos = 'title';
cfg.layout = '/lena13/home_users/users/hazem/fieldtrip-20161231/template/layout/neuromag306mag.lay'
figure; ft_topoplotTFR(cfg, power_low_contact);

% Topography beta
cfg = [];
cfg.xlim = [13 30] %[1 : 5 : 40];
cfg.zlim         = [-1.6e-28 1.6e-28];  %'maxmin' color for power
cfg.colormap =     jet ;
cfg.parameter    = 'powspctrm' 
%cfg.marker       = 'on';
cfg.comment = 'xlim';
cfg.commentpos = 'title';
cfg.layout = '/lena13/home_users/users/hazem/fieldtrip-20161231/template/layout/neuromag306mag.lay'
figure; ft_topoplotTFR(cfg, power_low_contact);

% Topography gamma
cfg = [];
cfg.xlim = [30 40] %[1 : 5 : 40];
cfg.zlim         = [-3.3e-29 3.3e-29];  %'maxmin' color for power
cfg.colormap =     jet ;
cfg.parameter    = 'powspctrm' 
%cfg.marker       = 'on';
cfg.comment = 'xlim';
cfg.commentpos = 'title';
cfg.layout = '/lena13/home_users/users/hazem/fieldtrip-20161231/template/layout/neuromag306mag.lay'
figure; ft_topoplotTFR(cfg, power_low_contact);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Power spectrum log10 on all the magnetometers averaged togethers

% cfg = []
% cfg.zlim = [-3e-27 3e-27];   %limite valeur Z (power)
% cfg.showlabels   = 'yes';
% cfg.layout = '/lena13/home_users/users/hazem/fieldtrip-20161231/template/layout/neuromag306mag.lay'   % for magnetometers only 
% %cfg.parameter     = field to be plotted on y-axis (default depends on data.dimord) 'avg', 'powspctrm' or 'cohspctrm'
% figure; ft_singleplotER(cfg, power_low_contact, power_high_contact);

clear all
input_filename = '/pclxserver/raid11/data/SELFI/Grand_average/High_contact_power_avgoverchannels.mat'
load (input_filename);
input_filename = '/pclxserver/raid11/data/SELFI/Grand_average/Low_contact_power_avgoverchannels.mat'
load (input_filename);

hold on
plot (power_high_contact_avgoverchannels.freq, log10(power_high_contact_avgoverchannels.powspctrm), 'r-')
plot (power_low_contact_avgoverchannels.freq, log10(power_low_contact_avgoverchannels.powspctrm), 'b-')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Individual Power Spectrum for all the sujects and each runs (not averaged over channels)

clear all

% load the data
load( '/pclxserver/raid11/data/SELFI/Grand_average/allsubj_power_High_contact_R1.mat');
load( '/pclxserver/raid11/data/SELFI/Grand_average/allsubj_power_High_contact_R2.mat');
load( '/pclxserver/raid11/data/SELFI/Grand_average/allsubj_power_Low_contact_R1.mat');
load( '/pclxserver/raid11/data/SELFI/Grand_average/allsubj_power_Low_contact_R2.mat');

% %Vizualisation not averaged over channels
% % freq = [5 15]; %Freq for the rectangle
% % ymax = log10(1e-26);
% % 
% % figure; 
% % for isub = 1:19
% %     subplot(4,5,isub) %Plot sur 4 lignes, 5 sujets par ligne,
% %     % use the rectangle to indicate the freq range of interest
% %     %rectangle('Position',[freq(1) 0 (freq(2)-freq(1)) ymax],'FaceColor',[0.7 0.7 0.7]); % [x y w h] in data units. The x and y elements determine the location and the w and h elements determine the size. The function plots into the current axes without clearing existing content from the axes.
% %     hold on;
% %     % plot the lines in front of the rectangle
% %     plot(allsubj_power_High_contact_R1{isub}.freq, log10(allsubj_power_High_contact_R1{isub}.powspctrm), 'r-')
% %     plot(allsubj_power_High_contact_R2{isub}.freq, log10(allsubj_power_High_contact_R2{isub}.powspctrm), 'r--')
% %     plot(allsubj_power_Low_contact_R1{isub}.freq, log10(allsubj_power_Low_contact_R1{isub}.powspctrm), 'b-')
% %     plot(allsubj_power_Low_contact_R2{isub}.freq, log10(allsubj_power_Low_contact_R2{isub}.powspctrm), 'b--')
% %     title(strcat('subject ',num2str(isub)))
% % %     ylim([0 log10(1e-26)]) 
% %     xlim([0 40])
% % end 
% % % subplot(5,5,21);   %Add the legend
% % % text(15,-26,'High_contact_run01','color','r') ;
% % % text(15,-27,'High_contact_run02 ---','color','r'); 
% % % text(15,-28,'Low_contact_run01','color','b') ;
% % % text(15,-29,'Low_contact_run02 ---','color','b')
% % % axis off
% 
% %Average over the subjects for each condition/run  
% 
% log_power = nan(102, 40)
% cfg = []
% cfg.parameter  = 'powspctrm'
% allsubj_power_High_contact_R1_avgoversubj = ft_freqgrandaverage(cfg, allsubj_power_High_contact_R1{:, 1})
% log_power = log10(allsubj_power_High_contact_R1_avgoversubj.powspctrm)
% allsubj_power_High_contact_R1_avgoversubj_log.powspctrm = allsubj_power_High_contact_R1_avgoversubj
% allsubj_power_High_contact_R1_avgoversubj_log.powspctrm = log_power
% allsubj_power_High_contact_R1_avgoversubj_log.label = allsubj_power_High_contact_R1_avgoversubj.label
% allsubj_power_High_contact_R1_avgoversubj_log.freq = allsubj_power_High_contact_R1_avgoversubj.freq 
% allsubj_power_High_contact_R1_avgoversubj_log.cfg = allsubj_power_High_contact_R1_avgoversubj.cfg
% allsubj_power_High_contact_R1_avgoversubj_log.dimord = allsubj_power_High_contact_R1_avgoversubj.dimord
% 
% log_power = nan(102, 40)
% cfg = []
% cfg.parameter  = 'powspctrm'
% allsubj_power_High_contact_R2_avgoversubj = ft_freqgrandaverage(cfg, allsubj_power_High_contact_R2{:, 1})
% log_power = log10(allsubj_power_High_contact_R2_avgoversubj.powspctrm)
% allsubj_power_High_contact_R2_avgoversubj_log.powspctrm = allsubj_power_High_contact_R2_avgoversubj
% allsubj_power_High_contact_R2_avgoversubj_log.powspctrm = log_power
% allsubj_power_High_contact_R2_avgoversubj_log.label = allsubj_power_High_contact_R2_avgoversubj.label
% allsubj_power_High_contact_R2_avgoversubj_log.freq = allsubj_power_High_contact_R2_avgoversubj.freq 
% allsubj_power_High_contact_R2_avgoversubj_log.cfg = allsubj_power_High_contact_R2_avgoversubj.cfg
% allsubj_power_High_contact_R2_avgoversubj_log.dimord = allsubj_power_High_contact_R2_avgoversubj.dimord
% 
% log_power = nan(102, 40)
% cfg = []
% cfg.parameter  = 'powspctrm'
% allsubj_power_Low_contact_R1_avgoversubj = ft_freqgrandaverage(cfg, allsubj_power_Low_contact_R1{:, 1})
% log_power = log10(allsubj_power_Low_contact_R1_avgoversubj.powspctrm)
% allsubj_power_Low_contact_R1_avgoversubj_log.powspctrm = allsubj_power_Low_contact_R1_avgoversubj
% allsubj_power_Low_contact_R1_avgoversubj_log.powspctrm = log_power
% allsubj_power_Low_contact_R1_avgoversubj_log.label = allsubj_power_Low_contact_R1_avgoversubj.label
% allsubj_power_Low_contact_R1_avgoversubj_log.freq = allsubj_power_Low_contact_R1_avgoversubj.freq 
% allsubj_power_Low_contact_R1_avgoversubj_log.cfg = allsubj_power_Low_contact_R1_avgoversubj.cfg
% allsubj_power_Low_contact_R1_avgoversubj_log.dimord = allsubj_power_Low_contact_R1_avgoversubj.dimord
% 
% log_power = nan(102, 40)
% cfg = []
% cfg.parameter  = 'powspctrm'
% allsubj_power_Low_contact_R2_avgoversubj = ft_freqgrandaverage(cfg, allsubj_power_Low_contact_R2{:, 1})
% log_power = log10(allsubj_power_Low_contact_R2_avgoversubj.powspctrm)
% allsubj_power_Low_contact_R2_avgoversubj_log.powspctrm = allsubj_power_Low_contact_R2_avgoversubj
% allsubj_power_Low_contact_R2_avgoversubj_log.powspctrm = log_power
% allsubj_power_Low_contact_R2_avgoversubj_log.label = allsubj_power_Low_contact_R2_avgoversubj.label
% allsubj_power_Low_contact_R2_avgoversubj_log.freq = allsubj_power_Low_contact_R2_avgoversubj.freq 
% allsubj_power_Low_contact_R2_avgoversubj_log.cfg = allsubj_power_Low_contact_R2_avgoversubj.cfg
% allsubj_power_Low_contact_R2_avgoversubj_log.dimord = allsubj_power_Low_contact_R2_avgoversubj.dimord
% 
% 
% %Multiplot_visu
% cfg = []
% cfg.channels = 'Meg*1'
% % cfg.zlim = [log(-3e-27) log(3e-27)];   %limite valeur Z (power)
% cfg.showlabels   = 'yes';
% cfg.layout = '/lena13/home_users/users/hazem/fieldtrip-20161231/template/layout/neuromag306mag.lay'   % for magnetometers only 
% cfg.parameter     = 'powspctrm'
% figure; ft_multiplotER(cfg, allsubj_power_Low_contact_R1_avgoversubj_log, allsubj_power_Low_contact_R2_avgoversubj_log, allsubj_power_High_contact_R1_avgoversubj_log, allsubj_power_High_contact_R2_avgoversubj_log);
% axis off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Individual Power Spectrum for all the subjects: COND*RUN (averaged over channels)

clear all

% load the data
load( '/pclxserver/raid11/data/SELFI/Grand_average/allsubj_power_High_contact_R1_avgoverchannels.mat');
load( '/pclxserver/raid11/data/SELFI/Grand_average/allsubj_power_High_contact_R2_avgoverchannels.mat');
load( '/pclxserver/raid11/data/SELFI/Grand_average/allsubj_power_Low_contact_R1_avgoverchannels.mat');
load( '/pclxserver/raid11/data/SELFI/Grand_average/allsubj_power_Low_contact_R2_avgoverchannels.mat');

%Vizualisation
freq = [5 15]; %Freq for the rectangle
ymax = log10(1e-26);
index_suj = [1, 2, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 19, 20, 21, 22]; 

figure; 
for isub = 1:19
    subplot(4,5,isub) %Plot sur 4 lignes, 5 sujets par ligne,
    % use the rectangle to indicate the freq range of interest
    %rectangle('Position',[freq(1) 0 (freq(2)-freq(1)) ymax],'FaceColor',[0.7 0.7 0.7]); % [x y w h] in data units. The x and y elements determine the location and the w and h elements determine the size. The function plots into the current axes without clearing existing content from the axes.
    hold on;
    % plot the lines in front of the rectangle
    plot(allsubj_power_High_contact_R1_avgoverchannels{isub}.freq, log10(allsubj_power_High_contact_R1_avgoverchannels{isub}.powspctrm), 'r-')
    plot(allsubj_power_High_contact_R2_avgoverchannels{isub}.freq, log10(allsubj_power_High_contact_R2_avgoverchannels{isub}.powspctrm), 'r--')
    plot(allsubj_power_Low_contact_R1_avgoverchannels{isub}.freq, log10(allsubj_power_Low_contact_R1_avgoverchannels{isub}.powspctrm), 'b-')
    plot(allsubj_power_Low_contact_R2_avgoverchannels{isub}.freq, log10(allsubj_power_Low_contact_R2_avgoverchannels{isub}.powspctrm), 'b--')
    title(strcat('subject ',num2str(index_suj(isub))))
%     ylim([0 log10(1e-26)]) 
    xlim([0 40])
end 
% subplot(5,5,21);   %Add the legend
% text(15,-26,'High_contact_run01','color','r') ;
% text(15,-27,'High_contact_run02 ---','color','r'); 
% text(15,-28,'Low_contact_run01','color','b') ;
% text(15,-29,'Low_contact_run02 ---','color','b')
% axis off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Individual Power Spectrum for all the subjects: COND (averaged over runs and channels)

clear all

% load the data
load( '/pclxserver/raid11/data/SELFI/Grand_average/allsubj_power_High_contact_avgoverchannels.mat');
load( '/pclxserver/raid11/data/SELFI/Grand_average/allsubj_power_Low_contact_avgoverchannels.mat');

%Vizualisation
% freq = [5 15]; %Freq for the rectangle
% ymax = log10(1e-26);

index_suj = [1, 2, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 19, 20, 21, 22]; 

figure; 
for isub = 1:19
    subplot(4,5,isub) %Plot sur 4 lignes, 5 sujets par ligne,
    % use the rectangle to indicate the freq range of interest
    %rectangle('Position',[freq(1) 0 (freq(2)-freq(1)) ymax],'FaceColor',[0.7 0.7 0.7]); % [x y w h] in data units. The x and y elements determine the location and the w and h elements determine the size. The function plots into the current axes without clearing existing content from the axes.
    hold on;
    % plot the lines in front of the rectangle
    plot(allsubj_power_High_contact_avgoverchannels{isub}.freq, log10(allsubj_power_High_contact_avgoverchannels{isub}.powspctrm), 'r-')
    plot(allsubj_power_Low_contact_avgoverchannels{isub}.freq, log10(allsubj_power_Low_contact_avgoverchannels{isub}.powspctrm), 'b-')
    
    title(strcat('subject ',num2str(index_suj(isub))))
%     ylim([0 log10(1e-26)]) 
    xlim([0 40])
end 
% subplot(5,5,21);   %Add the legend
% text(15,-26,'High_contact_run01','color','r') ;
% text(15,-27,'High_contact_run02 ---','color','r'); 
% text(15,-28,'Low_contact_run01','color','b') ;
% text(15,-29,'Low_contact_run02 ---','color','b')
% axis off

%Average over the subjects: H  
ligne = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38];

cfg = []
cfg.parameter  = 'powspctrm'
allsubj_power_High_contact_avgoverchannels_avgoversubj_avgoverruns = ft_freqgrandaverage(cfg, allsubj_power_High_contact_avgoverchannels{ligne, 1})

cfg = []
cfg.parameter  = 'powspctrm'
allsubj_power_Low_contact_avgoverchannels_avgoversubj_avgoverruns = ft_freqgrandaverage(cfg, allsubj_power_Low_contact_avgoverchannels{ligne, 1})


%Visu_mean_over_subj
figure;
hold on;
plot(allsubj_power_High_contact_avgoverchannels_avgoversubj_avgoverruns.freq, log10(allsubj_power_High_contact_avgoverchannels_avgoversubj_avgoverruns.powspctrm), 'r-')
plot(allsubj_power_Low_contact_avgoverchannels_avgoversubj_avgoverruns.freq, log10(allsubj_power_Low_contact_avgoverchannels_avgoversubj_avgoverruns.powspctrm), 'b-')
% ylim([0 log10(1e-26)]) 
xlim([0 40])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Vizu Power Spectrum on the subject with order HL and conditions taken separatly (averaged over channels and runs)

clear all

% load the data
load( '/pclxserver/raid11/data/SELFI/Grand_average/allsubj_power_High_contact_avgoverchannels.mat');
load( '/pclxserver/raid11/data/SELFI/Grand_average/allsubj_power_Low_contact_avgoverchannels.mat');

%Vizualisation
freq = [5 15]; %Freq for the rectangle
ymax = log10(1e-26);
ligne = [1, 2, 4, 5, 9, 10, 11, 16, 17];
num_graph = 1

figure; 
for isub = ligne
    subplot(2,5,num_graph) %Plot sur 2 lignes, 5 sujets par ligne,
    hold on;
    % plot the lines in front of the rectangle
    plot(allsubj_power_High_contact_avgoverchannels{isub}.freq, log10(allsubj_power_High_contact_avgoverchannels{isub}.powspctrm), 'r-')
    plot(allsubj_power_Low_contact_avgoverchannels{isub}.freq, log10(allsubj_power_Low_contact_avgoverchannels{isub}.powspctrm), 'b-')
%     ylim([0 log10(1e-26)]) 
    xlim([0 40])
    num_graph = num_graph + 1
end 

%Average over the subjects  
% ligne = [1, 2, 4, 5, 9, 10, 11, 16, 17];

cfg = []
cfg.parameter  = 'powspctrm'
allsubj_power_High_contact_avgoverchannels_avgoversubj = ft_freqgrandaverage(cfg, allsubj_power_High_contact_avgoverchannels{ligne, 1})

cfg = []
cfg.parameter  = 'powspctrm'
allsubj_power_Low_contact_avgoverchannels_avgoversubj = ft_freqgrandaverage(cfg, allsubj_power_Low_contact_avgoverchannels{ligne, 1})


%Visu_mean_over_subj
figure;
hold on;
plot(allsubj_power_High_contact_avgoverchannels_avgoversubj.freq, log10(allsubj_power_High_contact_avgoverchannels_avgoversubj.powspctrm), 'r-')
plot(allsubj_power_Low_contact_avgoverchannels_avgoversubj.freq, log10(allsubj_power_Low_contact_avgoverchannels_avgoversubj.powspctrm), 'b-')
% ylim([0 log10(1e-26)]) 
xlim([0 40])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Vizu Power Spectrum on the subject with order LH and conditions taken separatly (averaged over channels and runs)

clear all

% load the data
load( '/pclxserver/raid11/data/SELFI/Grand_average/allsubj_power_High_contact_avgoverchannels.mat');
load( '/pclxserver/raid11/data/SELFI/Grand_average/allsubj_power_Low_contact_avgoverchannels.mat');

%Vizualisation
freq = [5 15]; %Freq for the rectangle
ymax = log10(1e-26);
% ligne = [3, 6, 7, 8, 12, 13, 14, 15, 18, 19];
num_graph = 1

figure; 
for isub = ligne
    subplot(2,5,num_graph) %Plot sur 2 lignes, 5 sujets par ligne,
    hold on;
    % plot the lines in front of the rectangle
    plot(allsubj_power_High_contact_avgoverchannels{isub}.freq, log10(allsubj_power_High_contact_avgoverchannels{isub}.powspctrm), 'r-')
    plot(allsubj_power_Low_contact_avgoverchannels{isub}.freq, log10(allsubj_power_Low_contact_avgoverchannels{isub}.powspctrm), 'b-')
%     ylim([0 log10(1e-26)]) 
    xlim([0 40])
    num_graph = num_graph + 1
end 

%Average over the subjects  
% ligne =  [3, 6, 7, 8, 12, 13, 14, 15, 18, 19];

cfg = []
cfg.parameter  = 'powspctrm'
allsubj_power_High_contact_avgoverchannels_avgoversubj = ft_freqgrandaverage(cfg, allsubj_power_High_contact_avgoverchannels{ligne, 1})

cfg = []
cfg.parameter  = 'powspctrm'
allsubj_power_Low_contact_avgoverchannels_avgoversubj = ft_freqgrandaverage(cfg, allsubj_power_Low_contact_avgoverchannels{ligne, 1})


%Visu_mean_over_subj
figure;
hold on;
plot(allsubj_power_High_contact_avgoverchannels_avgoversubj.freq, log10(allsubj_power_High_contact_avgoverchannels_avgoversubj.powspctrm), 'r-')
plot(allsubj_power_Low_contact_avgoverchannels_avgoversubj.freq, log10(allsubj_power_Low_contact_avgoverchannels_avgoversubj.powspctrm), 'b-')
% ylim([0 log10(1e-26)]) 
xlim([0 40])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Vizu Power Spectrum on the subject with order HL and conditions/runs taken separatly (averaged over channels)

clear all

% load the data
load( '/pclxserver/raid11/data/SELFI/Grand_average/allsubj_power_High_contact_R1_avgoverchannels.mat');
load( '/pclxserver/raid11/data/SELFI/Grand_average/allsubj_power_High_contact_R2_avgoverchannels.mat');
load( '/pclxserver/raid11/data/SELFI/Grand_average/allsubj_power_Low_contact_R1_avgoverchannels.mat');
load( '/pclxserver/raid11/data/SELFI/Grand_average/allsubj_power_Low_contact_R2_avgoverchannels.mat');

%Vizualisation
freq = [5 15]; %Freq for the rectangle
ymax = log10(1e-26);
ligne = [1, 2, 4, 5, 9, 10, 11, 16, 17];
num_graph = 1

figure; 
for isub = ligne
    subplot(2,5,num_graph) %Plot sur 4 lignes, 5 sujets par ligne,
    % use the rectangle to indicate the freq range of interest
    %rectangle('Position',[freq(1) 0 (freq(2)-freq(1)) ymax],'FaceColor',[0.7 0.7 0.7]); % [x y w h] in data units. The x and y elements determine the location and the w and h elements determine the size. The function plots into the current axes without clearing existing content from the axes.
    hold on;
    % plot the lines in front of the rectangle
    plot(allsubj_power_High_contact_R1_avgoverchannels{isub}.freq, log10(allsubj_power_High_contact_R1_avgoverchannels{isub}.powspctrm), 'r-')
    plot(allsubj_power_High_contact_R2_avgoverchannels{isub}.freq, log10(allsubj_power_High_contact_R2_avgoverchannels{isub}.powspctrm), 'r--')
    plot(allsubj_power_Low_contact_R1_avgoverchannels{isub}.freq, log10(allsubj_power_Low_contact_R1_avgoverchannels{isub}.powspctrm), 'b-')
    plot(allsubj_power_Low_contact_R2_avgoverchannels{isub}.freq, log10(allsubj_power_Low_contact_R2_avgoverchannels{isub}.powspctrm), 'b--')
    title(strcat('subject ',num2str(isub)))
%     ylim([0 log10(1e-26)]) 
    xlim([0 40])
    num_graph = num_graph + 1
end 

%Average over the subjects  
ligne = [1, 2, 4, 5, 9, 10, 11, 16, 17];

cfg = []
cfg.parameter  = 'powspctrm'
allsubj_power_High_contact_R1_avgoverchannels_avgoversubj = ft_freqgrandaverage(cfg, allsubj_power_High_contact_R1_avgoverchannels{ligne, 1})

cfg = []
cfg.parameter  = 'powspctrm'
allsubj_power_High_contact_R2_avgoverchannels_avgoversubj = ft_freqgrandaverage(cfg, allsubj_power_High_contact_R2_avgoverchannels{ligne, 1})

cfg = []
cfg.parameter  = 'powspctrm'
allsubj_power_Low_contact_R1_avgoverchannels_avgoversubj = ft_freqgrandaverage(cfg, allsubj_power_Low_contact_R1_avgoverchannels{ligne, 1})

cfg = []
cfg.parameter  = 'powspctrm'
allsubj_power_Low_contact_R2_avgoverchannels_avgoversubj = ft_freqgrandaverage(cfg, allsubj_power_Low_contact_R2_avgoverchannels{ligne, 1})


%Visu_mean_over_subj
figure;
hold on;
plot(allsubj_power_High_contact_R1_avgoverchannels_avgoversubj.freq, log10(allsubj_power_High_contact_R1_avgoverchannels_avgoversubj.powspctrm), 'r-')
plot(allsubj_power_High_contact_R2_avgoverchannels_avgoversubj.freq, log10(allsubj_power_High_contact_R2_avgoverchannels_avgoversubj.powspctrm), 'r--')
plot(allsubj_power_Low_contact_R1_avgoverchannels_avgoversubj.freq, log10(allsubj_power_Low_contact_R1_avgoverchannels_avgoversubj.powspctrm), 'b-')
plot(allsubj_power_Low_contact_R2_avgoverchannels_avgoversubj.freq, log10(allsubj_power_Low_contact_R2_avgoverchannels_avgoversubj.powspctrm), 'b--')
% ylim([0 log10(1e-26)]) 
xlim([0 40])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Vizu Power Spectrum on the subject with order LH and conditions/runs taken separatly (averaged over channels)

clear all

% load the data
load( '/pclxserver/raid11/data/SELFI/Grand_average/allsubj_power_High_contact_R1_avgoverchannels.mat');
load( '/pclxserver/raid11/data/SELFI/Grand_average/allsubj_power_High_contact_R2_avgoverchannels.mat');
load( '/pclxserver/raid11/data/SELFI/Grand_average/allsubj_power_Low_contact_R1_avgoverchannels.mat');
load( '/pclxserver/raid11/data/SELFI/Grand_average/allsubj_power_Low_contact_R2_avgoverchannels.mat');

%Vizualisation
freq = [5 15]; %Freq for the rectangle
ymax = log10(1e-26);
ligne = [3, 6, 7, 8, 12, 13, 14, 15, 18, 19];
num_graph = 1

figure; 
for isub = ligne
    subplot(2,5,num_graph) %Plot sur 4 lignes, 5 sujets par ligne,
    % use the rectangle to indicate the freq range of interest
    %rectangle('Position',[freq(1) 0 (freq(2)-freq(1)) ymax],'FaceColor',[0.7 0.7 0.7]); % [x y w h] in data units. The x and y elements determine the location and the w and h elements determine the size. The function plots into the current axes without clearing existing content from the axes.
    hold on;
    % plot the lines in front of the rectangle
    plot(allsubj_power_High_contact_R1_avgoverchannels{isub}.freq, log10(allsubj_power_High_contact_R1_avgoverchannels{isub}.powspctrm), 'r-')
    plot(allsubj_power_High_contact_R2_avgoverchannels{isub}.freq, log10(allsubj_power_High_contact_R2_avgoverchannels{isub}.powspctrm), 'r--')
    plot(allsubj_power_Low_contact_R1_avgoverchannels{isub}.freq, log10(allsubj_power_Low_contact_R1_avgoverchannels{isub}.powspctrm), 'b-')
    plot(allsubj_power_Low_contact_R2_avgoverchannels{isub}.freq, log10(allsubj_power_Low_contact_R2_avgoverchannels{isub}.powspctrm), 'b--')
    title(strcat('subject ',num2str(isub)))
%     ylim([0 log10(1e-26)]) 
    xlim([0 40])
    num_graph = num_graph + 1
end 

%Average over the subjects  
ligne = [3, 6, 7, 8, 12, 13, 14, 15, 18, 19];

cfg = []
cfg.parameter  = 'powspctrm'
allsubj_power_High_contact_R1_avgoverchannels_avgoversubj = ft_freqgrandaverage(cfg, allsubj_power_High_contact_R1_avgoverchannels{ligne, 1})

cfg = []
cfg.parameter  = 'powspctrm'
allsubj_power_High_contact_R2_avgoverchannels_avgoversubj = ft_freqgrandaverage(cfg, allsubj_power_High_contact_R2_avgoverchannels{ligne, 1})

cfg = []
cfg.parameter  = 'powspctrm'
allsubj_power_Low_contact_R1_avgoverchannels_avgoversubj = ft_freqgrandaverage(cfg, allsubj_power_Low_contact_R1_avgoverchannels{ligne, 1})

cfg = []
cfg.parameter  = 'powspctrm'
allsubj_power_Low_contact_R2_avgoverchannels_avgoversubj = ft_freqgrandaverage(cfg, allsubj_power_Low_contact_R2_avgoverchannels{ligne, 1})


%Visu_mean_over_subj
figure;
hold on;
plot(allsubj_power_High_contact_R1_avgoverchannels_avgoversubj.freq, log10(allsubj_power_High_contact_R1_avgoverchannels_avgoversubj.powspctrm), 'r-')
plot(allsubj_power_High_contact_R2_avgoverchannels_avgoversubj.freq, log10(allsubj_power_High_contact_R2_avgoverchannels_avgoversubj.powspctrm), 'r--')
plot(allsubj_power_Low_contact_R1_avgoverchannels_avgoversubj.freq, log10(allsubj_power_Low_contact_R1_avgoverchannels_avgoversubj.powspctrm), 'b-')
plot(allsubj_power_Low_contact_R2_avgoverchannels_avgoversubj.freq, log10(allsubj_power_Low_contact_R2_avgoverchannels_avgoversubj.powspctrm), 'b--')
% ylim([0 log10(1e-26)]) 
xlim([0 40])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Power Spectrum Cond*Order (average over subjects)

clear all

% load the data
load( '/pclxserver/raid11/data/SELFI/Grand_average/allsubj_power_High_contact_avgoverchannels.mat');
load( '/pclxserver/raid11/data/SELFI/Grand_average/allsubj_power_Low_contact_avgoverchannels.mat');



%Vizualisation
% freq = [5 15]; %Freq for the rectangle
% ymax = log10(1e-26);

index_suj = [1, 2, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 19, 20, 21, 22]; 


%Average over the HL subjects:
ligne = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38];

cfg = []
cfg.parameter  = 'powspctrm'
allsubj_power_High_contact_avgoverchannels_avgoversubj_avgoverruns = ft_freqgrandaverage(cfg, allsubj_power_High_contact_avgoverchannels{ligne, 1})

cfg = []
cfg.parameter  = 'powspctrm'
allsubj_power_Low_contact_avgoverchannels_avgoversubj_avgoverruns = ft_freqgrandaverage(cfg, allsubj_power_Low_contact_avgoverchannels{ligne, 1})



load( '/pclxserver/raid11/data/SELFI/Grand_average/allsubj_power_High_contact_avgoverchannels.mat');
load( '/pclxserver/raid11/data/SELFI/Grand_average/allsubj_power_Low_contact_avgoverchannels.mat');

%Average over the LH subjects:
ligne =  [3, 6, 7, 8, 12, 13, 14, 15, 18, 19];

cfg = []
cfg.parameter  = 'powspctrm'
allsubj_power_High_contact_avgoverchannels_avgoversubj = ft_freqgrandaverage(cfg, allsubj_power_High_contact_avgoverchannels{ligne, 1})

cfg = []
cfg.parameter  = 'powspctrm'
allsubj_power_Low_contact_avgoverchannels_avgoversubj = ft_freqgrandaverage(cfg, allsubj_power_Low_contact_avgoverchannels{ligne, 1})


%Visu_mean_over_subj
figure;
hold on;
plot(allsubj_power_High_contact_avgoverchannels_avgoversubj_avgoverruns.freq, log10(allsubj_power_High_contact_avgoverchannels_avgoversubj_avgoverruns.powspctrm), 'r-')
plot(allsubj_power_Low_contact_avgoverchannels_avgoversubj_avgoverruns.freq, log10(allsubj_power_Low_contact_avgoverchannels_avgoversubj_avgoverruns.powspctrm), 'b-')
% ylim([0 log10(1e-26)]) 
xlim([0 40])