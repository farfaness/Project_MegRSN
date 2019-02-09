%% To display informations on the header of the filename

hdr = ft_read_header(filename); % change the filename
disp(hdr);
disp(hdr.label');
disp(hdr.chantype');
disp(hdr.chanunit');

%% To visualize the data

cfg = [];
cfg.dataset = filename;
cfg.channel = 'MEG';    % or 'Megmag' or 'Meggrad'
cfg.viewmode = 'vertical';
ft_databrowser(cfg);