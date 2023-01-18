% spike-sorting script
function spikeSort_maze(pt,exp)

DEBUG = 0;

dropboxLink()

header = getmemMazeExperimentHeader(pt,exp);
datasetFilename = fullfile(header.processedDataPath,sprintf('%s_EXP%d_dataset.mat',header.id,header.experimentNum));
a = dir(datasetFilename);
if ~isempty(a)
    
    load(datasetFilename)
    load(fullfile(header.montagePath))
end

whatToDo = questdlg('Select input source?','source','AveragedRef','RAW','RAW');
%% Spike sorting on averaged channels  -----------------------------------------------------------------
if isfield(header,'NS5ChannelLabels') % Spikes recorded on BR
    channels = header.NS5ChannelLabels;
    channels(EXP_DATA.ELECTRODES*8+1:end) = [];
    header.spikes_sampling_Rate_Hz = header.MICRO_Sampling_Rate_Hz;
else
    channels = 1:EXP_DATA.ELECTRODES*8;
end
if strcmp(whatToDo,'RAW')
    source_dir = header.spikesDataPath;
else
    source_dir = header.processed_AverageRef;
end
RERUN = 0;
PRINT_ONLY = 1;
cd(source_dir)

for jj = 1:EXP_DATA.ELECTRODES
    channels = (jj-1)*8+1:(jj)*8;
    indices =  1:length(channels);
    for ii = indices
        
        disp(sprintf('working on channel %d',channels(ii)))
        
        filename = fullfile(source_dir,sprintf('CSC%d.mat',channels(ii)));
        filename2 = fullfile(source_dir,sprintf('CSC%d_001.mat',channels(ii)));
        if (isempty(dir(filename)) && isempty(dir(filename2)))
            disp(sprintf('%s doesn''t exists, skipped',filename))
            continue
        end
        if ~RERUN
            filename = fullfile(source_dir,sprintf('CSC%d_spikes.mat', channels(ii)));
            if ~isempty(dir(filename))
                disp(sprintf('%s already exists, skipped',filename))
                continue
            end
        end
        if ~isempty(dir(filename2))
            mergeFiles = 1;
            fprintf('using merged spike-files for sorting')
        else
            mergeFiles = 0;
        end
        
        detectionType = 'pos';
        getSpikesMaya(source_dir, channels(ii),mergeFiles, header.spikes_sampling_Rate_Hz,detectionType)
                                          % Notice SR should be updated for NLX files
        close all
        
    end
    
    EDIT_SPIKE_FILES = 1;
    spikeCoincidenceDetection(source_dir, channels, EDIT_SPIKE_FILES)
    
    header.SpikeSortingParams = []; % Need to load these values
    
    
    for ii = indices
        
        disp(sprintf('working on channel %d',channels(ii)))
        
        filename = fullfile(source_dir,sprintf('CSC%d.mat',channels(ii)));
        filename2 = fullfile(source_dir,sprintf('CSC%d_001.mat',channels(ii)));
        if (isempty(dir(filename)) && isempty(dir(filename2)))
            disp(sprintf('%s doesn''t exists, skipped',filename))
            continue
        end
        if ~RERUN
            filename = fullfile(source_dir,sprintf('times_CSC%d.mat', channels(ii)));
            if ~isempty(dir(filename))
                disp(sprintf('%s already exists, skipped',filename))
                if PRINT_ONLY
                    % TBD - code not working
             %       printClusteringResults(source_dir, channels(ii), header.spikes_sampling_Rate_Hz)
                end
                continue
            end
        end
        if ~isempty(dir(filename2))
            mergeFiles = 1;
            fprintf('using merged spike-files for sorting')
        else
            mergeFiles = 0;
        end
        
        doClusteringMaya(source_dir, channels(ii),header.spikes_sampling_Rate_Hz); % Notice SR should be updated for NLX files
        close all
        
    end
    
end

% prepare for manual sorting
for ii = 1:length(channels)
    ch = channels(ii);
    f1 = sprintf('CSC%d_001.mat',ch);
    f2 = sprintf('CSC%d.mat',ch);
    disp(sprintf('copy %s to %s and remove NaNs',f1,f2))
    copyfile(f1,f2)
    M1 = matfile(f2,'Writable',true);
    dataTemp = M1.data;
    indices = find(isnan(dataTemp));
    dataTemp(indices) = 0;
    M1.data = dataTemp;
    clear M1 dataTemp f1 f2
end

if DEBUG
%% Compare original data to avergedREf
figure
MICRO_SR = 32e3;
buffer = 10*MICRO_SR;
shift = 100;
% pt 496
rawData = 'H:\data_SLP5_backup\p496\EXP8\spikeSorting';
averagedRef = 'C:\Users\Maya\Data_p\p496\EXP8\averagedRef\';

% pt 488
averagedRef = 'J:\Data_p\p488\EXP4\averagedRef';
times_files_link =  fullfile(averagedRef,'spikeSortingResults_waveClus');
rawData =   'H:\Data_p\p488\EXP4\spikeSorting';
times_files_link_real = rawData;

% prepare for manual sorting
for ii = 1:length(channels)
    ch = channels(ii);
    f1 = fullfile(rawData,sprintf('CSC%d_001.mat',ch));
    f2 = fullfile(averagedRef,sprintf('CSC%d_001.mat',ch));
    f3 = fullfile(times_files_link_avg,sprintf('times_CSC%d.mat',ch));
    f4 = fullfile(times_files_link_real,sprintf('times_CSC%d.mat',ch));
    
    M1 = matfile(f1);
    M2 = matfile(f2);
    M3_sp = matfile(f3);
    M4_sp = matfile(f4);
    indices_sp = find( M4_sp.cluster_class(:,1) == 1);
    indices_sp_t = M4_sp.cluster_class(indices_sp(1),2);
    indices_sp_t_MICRO_SR = floor(header.MICRO_Sampling_Rate_Hz * indices_sp_t/1000);
    
    t_min = (1/(60*MICRO_SR))*(1:buffer);
    t0 = max(1,indices_sp_t_MICRO_SR-buffer);
    plot(M2.data(1,t0:indices_sp_t_MICRO_SR+buffer)+shift*(ii-1),'r')
    hold all
    line(ones(1,2)*indices_sp_t_MICRO_SR,get(gca,'ylim'),'color','r')
    plot(M1.data(1,t0:indices_sp_t_MICRO_SR+buffer)+shift*(ii-1),'b')
    hold all
    plot(M2.REF_average(1,t0:indices_sp_t_MICRO_SR+buffer)+shift*(ii-1)+200,'k')
    title('red - data-ref, blue - raw data, black - avg REF')
    
    hold all
    clear M1 
end
end

disp('now run manual wave_clus and update XLS')

% Delete merged files
% for ii = 1:length(channels)
%     ch = channels(ii);
%     f2 = sprintf('CSC%02d.mat',ch);
%     disp(sprintf('deleting %s',f2))
%     delete(f2)
% end

whatToDo = questdlg('Would you like to update spike-sorted results now?','SpikeSort','YES','NO','NO');
disp(sprintf('spike sorted list will be extracted from %s\\%s',header.excel_sheet,header.spikesorting_sheet))
disp(sprintf('spike sorted files will be extracted from %s',source_dir))
if strcmp(whatToDo,'YES')
    List = getSpikeSortedCellList(header.excel_sheet,header.spikesorting_sheet,header.N_spikeSortedCells);
    
    spikeSortingDir = source_dir;
    
    whatToDo = questdlg('remove artifacts appearing on all electrodes','artifacts','YES','NO','NO');
    if strcmp(whatToDo, 'YES')
        EDIT_SPIKE_FILES = 1;
        for ii = 1:EXP_DATA.ELECTRODES
            channels = (ii-1)*8+1:(ii-1)*8+8;
            spikeCoincidenceDetection(spikeSortingDir, channels, EDIT_SPIKE_FILES)
        end
    end
    
    % Load timestamps of spike timing for all cells
    clear spike_timestamps micro_channels_spike_summary
    spike_timestamps = []; spike_shapes_mean = []; spike_shapes_std = [];
    for ii = 1:length(List)
        channel = List(ii).Channel;
        cluster = List(ii).Cluster;
        filename = fullfile(spikeSortingDir, sprintf('times_CSC%d.mat',channel));
        cluster_class = load(filename,'cluster_class');
        cluster_class = cluster_class.cluster_class;
        spike_timestamps{ii} = cluster_class(cluster_class(:,1) == cluster,2); % msec
        spikeShapes_all = load(filename,'spikes');spikeShapes_all = spikeShapes_all.spikes;
        spikeShapes_SU = spikeShapes_all(cluster_class(:,1) == cluster,:);
        spike_shapes_mean{ii} = mean(spikeShapes_SU);
        spike_shapes_std{ii} = std(spikeShapes_SU);
    end
    
    micro_channels_spike_summary.spikesDataPath = spikeSortingDir;
    micro_channels_spike_summary.spike_timestamps = spike_timestamps;
    micro_channels_spike_summary.unit_list_xls = List;
    micro_channels_spike_summary.xls_filename = header.excel_sheet;
    micro_channels_spike_summary.spikesorting_sheet = header.spikesorting_sheet;
    micro_channels_spike_summary.spike_shapes_mean = spike_shapes_mean;
    micro_channels_spike_summary.spike_shapes_std = spike_shapes_std;
        
    save(fullfile(spikeSortingDir,sprintf('%s_spike_timestamps_post_processing',header.id)),'micro_channels_spike_summary')
end
