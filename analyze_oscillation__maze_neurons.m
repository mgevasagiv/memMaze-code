% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Analyze data to find low-frequency oscilaltions (lower then the theta
% band) around 1Hz 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% =========================================================================
% Maya - CA1 - interneurons - flight (only LIGHT) - egyptian fruit bat
% =========================================================================

%%
% clear all;

data_dir = 'E:\MAZE\data_p\singleUnit';
cell_IX = [16   131   197   203   205   206   209   210   212   215   216   221   225   228   230   231   233   234   235   238   245   249   251   253   256   271   273   275];  % FR>5Hz & spike_width < 0.28

data = {};
for ii_cell = 1:length(cell_IX)
    
    cell_num = cell_IX(ii_cell);
    
    file_name = sprintf('density_per_dim_statistics_cell_%03d', cell_num );
    file_to_load = fullfile(data_dir, file_name );
    clear density_per_dim_statistics_one_cell;
    load(file_to_load);
    
    SPIKES_Timestamps_temp = [];
    if length(density_per_dim_statistics_one_cell) == 3
        sessions_to_analyze = [1 3];
    else
        sessions_to_analyze = 1:length(density_per_dim_statistics_one_cell);
    end
    for session = sessions_to_analyze
        if isempty(density_per_dim_statistics_one_cell{session})
            continue;
        end
        SPIKES_Timestamps_temp = [  SPIKES_Timestamps_temp density_per_dim_statistics_one_cell{session}.PF_density_per_dim.spike_timestamps{3}];
    end

    data{ii_cell}.SPIKES_Timestamps = SPIKES_Timestamps_temp;
end
clear density_per_dim_statistics_one_cell;

res_dir = 'P:\Results\maya_CA1_interneurons_EFB\light\';
figure_dir = fullfile(res_dir, 'figures');
mkdir(figure_dir);
mkdir(res_dir);

%% define parameters
Time_around_temporal_autocorrelation_central_peak_in_msec = 2500;
bin_size_for_temporal_autocorrelogram = 10; % in msec
% bin_size_for_temporal_autocorrelogram = 50; % in msec
main_freqs = [1:1];
band_list = repmat(main_freqs, 2, 1)' + repmat([-0.5 0.5], length(main_freqs),1);
theta_band = [5 11];
zeroPadFactor = 16;
minimal_freq_val_for_reamining_spectrum_range = 0;
freq_around_theta_band_peak_in_Hz  = 0;
% freq_around_theta_band_peak_in_Hz  = 1; % As in Boccara Nature Neuroscience 2010 - Note that for a cell to be defined theta rhytmic
% % The ratio has to be > 5!!!
smoothing_window_in_Hz = 1;

shuffle_num_rep = 1;
spike_randomization_option = 1;
temporal_jitter_Sec = 10; % Time in Sec of jitter around spike occurance on EACH side
num_cells_to_analyze = Inf; % use this if you want just a partial analysis to check parameters on few cells...
run_shuffling = 0;

save_data = 1;

%% run the wrapper script
analyze_oscillation__wrapper







%%