function besa2edf_local(filename_in)

hdr = ft_read_header(filename_in);

cfg = [];
cfg.dataset = filename_in;
data = ft_preprocessing(cfg);

filename_out = [filename_in(1:end-4) '.edf'];
ft_write_data(filename_out, data.trial{1}, 'header', hdr, 'format', 'edf');
