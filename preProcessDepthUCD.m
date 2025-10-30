function preProcessDepthUCD(subj,exp)

dropbox_link = 'E:\Dropbox\';
% addpath('E:\Dropbox\Code\Nir Lab\Work\ONLINE_EEG_TRACKING_analysis')
% Load header + EXP data once it's ready
header = getmemMazeExperimentHeader(sprintf('%s',subj), exp);
% PtList = importXLSMazePatientList('E:\NeurDrift\NeurDriftSessions.xlsx'...
%             ,'Subjects',10);

PtList = readtable('E:\Dropbox\RanganathLab\work\mazeStudy\mazeSessions.xlsx');   % For tabular data

ptInd = [];

for ii = 1:length(PtList.subj)
    if strcmp(PtList.subj(ii),subj)
        if PtList.exp(ii) == exp
            ptInd = ii;
            break;
        end
    end
end
if isempty(ptInd)
    disp('pt entry not found in xls, searching general list')
    PtList = importXLSClosedLoopPatientList(fullfile(dropbox_link,'RanganathLab\work\mazeStudy\mazeSessions.xlsx')...
        ,'subjects',25);
    ptInd = [];
    for ii = 1:length(PtList)
        if strcmp(PtList(ii).subj,subj);
            if PtList(ii).exp == exp
                ptInd = ii;
                break;
            end
        end
    end
    PtEntry = PtList(ptInd);
    
else
    PtEntry = PtList(ptInd);
end
datasetFilename = fullfile(header.processedDataPath,sprintf('p%s_EXP%d_dataset.mat',header.id,header.experimentNum));
a = dir(datasetFilename);
if ~isempty(a)
    
    load(datasetFilename)
    
    % Load MacroMontage
    % load(fullfile(header.macroMontagePath,'MacroMontage.mat'))
    
    % Load spike data once it's available
    % load(fullfile(header.processed_AverageRef,sprintf('%s_spike_timestamps_post_processing',header.id)))
    
else % prepare dataset file
    
    if prep_macroMontage
        macroMontageSheet = 'MacroMontage';
        p_in.numeric_fields = {'Channel'}; % numeric fields in excel
        MacroMontage= [];
        
        if strfind(header.id,'p') % manual montage for UCLA patients
            cell_rows = 1:250;
            % MacroMontage = read_excel_sheet(header.excel_sheet,macroMontageSheet,cell_rows,p_in.numeric_fields,p_in);
            MacroMontage = readtable(header.excel_sheet);   % For tabular data

            % get rid of empty entries
            ii = length(MacroMontage.Channel);
            while(1)
                if isempty(MacroMontage.Area(ii))
                    MacroMontage(ii) = [];
                    disp(sprintf('removed entry %d',ii))
                else
                    break
                end
                ii = ii - 1;
            end
            if isempty(dir(header.root_data))
                mkdir(header.root_data)
            end
            save(header.macroMontagePath,'MacroMontage')
        end
        
    end
    %% This should be ran once for each session and saved for later usage
    %% Manually update the following few lines -
    %     header.excel_sheet = fullfile(dropbox_link,'Nir_Lab','Work', PtEntry.excel_sheet);
    %     header.spikesorting_sheet = PtEntry.spikesorting_sheet; % SpikeSorting Sheet in patient's xls
    %     header.N_spikeSortedCells = PtEntry.SpikeSortedUnits;
    %     if strcmp(PtEntry.spikesRecordingSys,'BR')
    %         header.montage_sheet = sprintf('MicroMontage_EXP%d',PtEntry.Nsessions);
    %     elseif strcmp(PtEntry.spikesRecordingSys,'NLX')
    %         header.montage_sheet = sprintf('MicroMontage',PtEntry.Nsessions);
    %     end
    if exp == 1
        timestampsTable = readtable(fullfile('E:\MAZE\ptData\',sprintf('%s_timestamp_table.csv',subj)));
    else
        timestampsTable = readtable(fullfile('E:\MAZE\ptData\',sprintf('%s_e%d_timestamp_table.csv',subj,exp)));
    end
    
    EXP_DATA.LINE_FREQUENCY = 60;
    
    for ii = 1:3
        rows = find(ismember(timestampsTable.Type , 'X') &ismember(timestampsTable.Repetition , ii) &...
            ~ismember(timestampsTable.Block , [10 11 12]));
        EXP_DATA.timestamps_us.X{ii} = timestampsTable.Time_onset(rows);
        
        rows = find(ismember(timestampsTable.Type , 'G') & ismember(timestampsTable.Repetition , ii)&...
            ~ismember(timestampsTable.Block , [10 11 12]));
        EXP_DATA.timestamps_us.G{ii} = timestampsTable.Time_onset(rows);
        
        rows = find(ismember(timestampsTable.Type , 'D') & ismember(timestampsTable.Repetition , ii)&...
            ~ismember(timestampsTable.Block , [10 11 12]));
        EXP_DATA.timestamps_us.D{ii} = timestampsTable.Time_onset(rows);
        
        rows = find(ismember(timestampsTable.Type , 'N') & ismember(timestampsTable.Repetition , ii)&...
            ~ismember(timestampsTable.Block , [10 11 12]));
        EXP_DATA.timestamps_us.N{ii} = timestampsTable.Time_onset(rows);
    end
    
    for ii = 1:3
        rows = find(ismember(timestampsTable.Type , 'X') &ismember(timestampsTable.Repetition , ii) &...
            ismember(timestampsTable.Block , [10 11 12]));
        EXP_DATA.timestamps_us.X_B4{ii} = timestampsTable.Time_onset(rows);
        
        rows = find(ismember(timestampsTable.Type , 'G') & ismember(timestampsTable.Repetition , ii)&...
            ismember(timestampsTable.Block , [10 11 12]));
        EXP_DATA.timestamps_us.G_B4{ii} = timestampsTable.Time_onset(rows);
        
        rows = find(ismember(timestampsTable.Type , 'D') & ismember(timestampsTable.Repetition , ii)&...
            ismember(timestampsTable.Block , [10 11 12]));
        EXP_DATA.timestamps_us.D_B4{ii} = timestampsTable.Time_onset(rows);
        
        rows = find(ismember(timestampsTable.Type , 'N') & ismember(timestampsTable.Repetition , ii)&...
            ismember(timestampsTable.Block , [10 11 12]));
        EXP_DATA.timestamps_us.N_B4{ii} = timestampsTable.Time_onset(rows);
    end
    
    str = 'test';
    testing_indices = find(~isnan(timestampsTable.Test_contextual_success));
    correctTest = timestampsTable.Test_contextual_success(testing_indices) == 1;
    inCorrectTest = timestampsTable.Test_contextual_success(testing_indices) == 0;
    onset_testing_c = timestampsTable.Time_onset(testing_indices(correctTest)-1);
    offset_testing_c = timestampsTable.Time_offset(testing_indices(correctTest)-1);
    onset_set = onset_testing_c;
    offset_set = offset_testing_c;
    
    EXP_DATA.onset_testing_correct = onset_set;
    EXP_DATA.offset_testing_correct = offset_set;
    
    str = 'nc_test';
    NC_testing_indices = find(~isnan(timestampsTable.Test_non_contextual_success));
    correctTest = timestampsTable.Test_non_contextual_success(NC_testing_indices) == 1;
    onset_NC_testing_c = timestampsTable.Time_onset(NC_testing_indices(correctTest)-1);
    offset_NC_testing_c = timestampsTable.Time_offset(NC_testing_indices(correctTest)-1);
    onset_set = onset_NC_testing_c;
    offset_set = offset_NC_testing_c;
    EXP_DATA.onset_testing_incorrect = onset_set;
    EXP_DATA.offset_testing_incorrect = offset_set;

    % Focus on memorized mazes
    rows = find(ismember(timestampsTable.Test_contextual_success ,1) );
    EXP_DATA.succsesful_maze = unique(timestampsTable.Maze(rows));
    
    
    filename = fullfile(header.processedDataPath,sprintf('%s_EXP%d_dataset.mat',header.id,header.experimentNum));
    save(filename,'header','EXP_DATA')
    
end

return

whatToDo = questdlg('Would you like to extract data files?','extract?','NS3','NS5','NO','NO');
if strcmp(whatToDo,'NS3')
    
    SR = NSx3Header.MetaTags.SamplingFreq;
    channels = 1:NSx3Header.MetaTags.ChannelCount;
    channel_labels = NSx3Header.MetaTags.ChannelID;
    output_folder = fullfile(header.processedDataPath,'LP_channels_2000Hz');
    
    for ii = channels
        
        disp(sprintf('working on channel %d',channels(ii)))
        channel_id = NSx3Header.MetaTags.ChannelID(channels(ii));
        CreateDataFilesFromNsx_LA(nsxPath3, channels(ii), channel_id, output_folder, header.Digital_to_uV, SR)
        
    end
    
elseif strcmp(whatToDo,'NS5')
    nsxPath5 = [header.nsFilePath{Index.MICRO},'.',header.nsFormat{Index.MICRO}];
    channels = 1:length(NSx5Header.MetaTags.ChannelID);
    channels(channels > EXP_DATA.ELECTRODES*8) = [];
    
    for ii = 1:length(channels)
        
        disp(sprintf('working on channel %d',channels(ii)))
        channel_id = NSx5Header.MetaTags.ChannelID(channels(ii));
        CreateDataFilesFromNsx_LA(nsxPath5, channels(ii), channel_id, header.spikesDataPath, ...
            header.Digital_to_uV, header.MICRO_Sampling_Rate_Hz)
        close all
        
    end
    
    %% Generate denoised with an average reference for ripple detection (NO High pass filtering)
    whatToDo = questdlg('Would you like to calc averaged ref for all electrodes?','REF?','ALL','PARTIAL','ALL');
    
    source_folder = header.spikesDataPath;
    target_folder = header.processed_AverageRef;
    if (isempty(dir(target_folder)))
        mkdir(target_folder)
    end
    if strcmp(whatToDo,'ALL')
        if isfield(header,'NS5ChannelLabels')
            channels = header.NS5ChannelLabels;
            channels(EXP_DATA.ELECTRODES*8+1:end) = [];
        else
            channels = 1:EXP_DATA.ELECTRODES*8;
        end
    else
        answer = input('enter numbers - ')
    end
    preProcessMICRO_avgRef()
    
end

whatToDo = questdlg('Would you like to create 1kHz data-files?','downsample','YES','NO','NO');
if strcmp(whatToDo,'YES')
    
    %% Generate a resampled version, to allow easier loading to EDF file
    source_folder = header.spikesDataPath;
    target_folder = header.processed_MICRO;
    if (isempty(dir(target_folder)))
        mkdir(target_folder)
    end
    
    channels = header.NS5ChannelLabels;
    channels(EXP_DATA.ELECTRODES*8+1:end) = [];
    header.MICRO_Sampling_Rate_Hz = PtEntry.samplingFreq;
    
    if isfield(header,'MontageCell')
        microMontage = header.MontageCell;
    else
        mm = matfile(header.montagePath);
        microMontage = mm.Montage;
    end
    for ii_c = 1:length(channels)
        
        channel_id = channels(ii_c);
        if isempty(dir(fullfile(source_folder,sprintf('CSC%d.mat',channel_id) )))
            fprintf('channel %d skipped\n',channel_id)
            continue
        end
        mLink = matfile(fullfile(source_folder,sprintf('CSC%d.mat',channel_id) ));
        % data = resample(mLink.data,Consts.DOWNSAMPLED_FREQUENCY,header.MICRO_Sampling_Rate_Hz);
        data = accurateResampling(mLink.data,header.MICRO_Sampling_Rate_Hz,Consts.DOWNSAMPLED_FREQUENCY); % July 2018 - upgrading to better resampling
        
        filename = fullfile(target_folder,sprintf('CSC%d.mat',channel_id));
        LocalHeader.samplingRate = Consts.DOWNSAMPLED_FREQUENCY;
        LocalHeader.channel_id = channel_id;
        LocalHeader.label = microMontage(channel_id).Area;
        LocalHeader.samplerate = Consts.DOWNSAMPLED_FREQUENCY;
        save(filename,'data','header','LocalHeader')
    end
    
    source_folder = fullfile(header.processedDataPath,'LP_channels_2000Hz');
    channels = header.NS3ChannelLabels;
    samplingRate = 2e3;
    for ii_c = 1:length(channels)
        
        channel_id = channels(ii_c);
        if isempty(dir(fullfile(source_folder,sprintf('CSC%d.mat',channel_id) )))
            fprintf('channel %d skipped',channel_id)
            continue
        end
        load(fullfile(source_folder,sprintf('CSC%d.mat',channel_id) ))
        % resampled_data = resample(data,Consts.DOWNSAMPLED_FREQUENCY,samplingRate);
        data = accurateResampling(data,samplingRate,Consts.DOWNSAMPLED_FREQUENCY); % July 2018 - upgrading to better resampling
        data = resampled_data;
        filename = fullfile(target_folder,sprintf('CSC%d.mat',channel_id));
        LocalHeader.samplingRate = Consts.DOWNSAMPLED_FREQUENCY;
        LocalHeader.channel_id = channel_id;
        for jj = 1:length(header.MontageCell)
            if header.MontageCell{jj,1} == channels(ii_c)
                LocalHeader.label = header.MontageCell{jj,2};
                break
            end
        end
        LocalHeader.samplerate = Consts.DOWNSAMPLED_FREQUENCY;
        save(filename,'data','header','LocalHeader')
    end
end

whatToDo = questdlg('Would you like SpikeSort Now?','SpikeSort','YES','NO','NO');
if strcmp(whatToDo,'YES')
    whatToDo = questdlg('Select input source?','source','AveragedRef','RAW','RAW');
    %% Spike sorting on averaged channels  -----------------------------------------------------------------
    if isfield(header,'NS5ChannelLabels')
        channels = header.NS5ChannelLabels;
        channels(EXP_DATA.ELECTRODES*8+1:end) = [];
    else
        channels = 1:EXP_DATA.ELECTRODES*8;
    end
    if strcmp(whatToDo,'RAW')
        source_dir = header.spikesDataPath;
    else
        source_dir = header.processed_AverageRef;
    end
    RERUN = 1;
    
    header.MICRO_Sampling_Rate_Hz = PtEntry.spikeSamplingRate_Hz;
    indices =  1:length(channels);
    parfor ii = indices
        
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
                continue
            end
        end
        if ~isempty(dir(filename2))
            mergeFiles = 1;
            fprintf('using merged spike-files for sorting')
        else
            mergeFiles = 0;
        end
        
        getSpikesMaya(source_dir, channels(ii),mergeFiles,header.MICRO_Sampling_Rate_Hz ); % Notice SR should be updated for NLX files
        doClusteringMaya(source_dir, channels(ii),header.MICRO_Sampling_Rate_Hz ); % Notice SR should be updated for NLX files
        close all
        
    end
    
    EDIT_SPIKE_FILES = 1;
    spikeSortingDir = source_dir;
    spikeCoincidenceDetection(spikeSortingDir, channels, EDIT_SPIKE_FILES)
    
    header.SpikeSortingParams = []; % Need to load these values
    
end
disp('Update single-units in pt''s XLS')

whatToDo = questdlg('Would you like to update spike-sorted results now?','SpikeSort','YES','NO','NO');
if strcmp(whatToDo,'YES')
    List = getSpikeSortedCellList(header.excel_sheet,header.spikesorting_sheet,header.N_spikeSortedCells);
    
    whatToDo = questdlg('Choose spike sorting location','SpikeSort',header.processed_AverageRef,header.spikesDataPath,header.processed_AverageRef);
    spikeSortingDir = whatToDo;
    
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
    spike_timestamps = [];
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

whatToDo = questdlg('Would you like to generate EDF file now?','EDF','YES','NO','NO');
if strcmp(whatToDo,'YES')
    
    whatToDo = questdlg('Which signal set to use?','EDF','MICRO','MACRO','MICRO');
    if strcmp(whatToDo,'MICRO')
        Source_folder_EDF =  header.processed_MICRO;
        uMontage = load(header.montagePath,'Montage'); uMontage = uMontage.Montage;
        channel_vec = 1:length(uMontage);
        str = 'MICRO';
    else
        Source_folder_EDF =  header.processed_MACRO;
        MacroMontage = load(header.macroMontagePath,'MacroMontage'); MacroMontage = MacroMontage.MacroMontage;
        channel_vec = 1:length(MacroMontage);
        EEG_SIG = [{'C3'},{'C4'},{'Pz'},{'Ez'}];
        str = 'MACRO';
    end
    clear DATA_MAT_EDF EDFHeader
    
    for ii_c = 1:length(channel_vec)
        
        channel_n = channel_vec(ii_c);
        disp(sprintf('loading channel #%d/%d', ii_c, length(channel_vec)))
        
        clear data
        if strcmp(whatToDo,'MICRO')
            channel = header.NS5ChannelLabels(channel_n);
        else
            channel = MacroMontage(channel_n).Channel;
        end
        filename = fullfile(Source_folder_EDF, sprintf('CSC%d.mat',channel));
        
        if isempty(dir(filename))
            fprintf('channel %d missing\n',channel) % These were inputs to closed loop
            continue
        else
            data = load(filename,'data'); % 30kHz downsampled channels
            data = data.data;
        end
        if ii_c == 1
            L = length(data);
        else
            if length(data) ~= L
                disp('%s length is uncompatible - dropped')
                data = zeros(L,1);
            end
        end
        
        DATA_MAT_EDF(ii_c,:) = data(:)';
        if strcmp(whatToDo,'MICRO')
            EDFHeader.labels{ii_c} = sprintf('%d-%s',channel,header.MontageCell{channel_n,2});
        else
            EDFHeader.labels{ii_c} = sprintf('MACRO-%d-%s',channel,MacroMontage(channel_n).Area);
        end
    end
    
    % Add MACRO channels - bypassed from clinical system
    channel_vec = 1:length(header.NS3ChannelLabels);
    for ii_cc = 1:length(channel_vec)
        
        clear data
        channel = header.NS3ChannelLabels(ii_cc);
        filename = fullfile(Source_folder_EDF, sprintf('CSC%d.mat',channel));
        
        if isempty(dir(filename))
            fprintf('channel %d missing\n',channel) % These were inputs to closed loop
            continue
        else
            data = load(filename,'data'); % 30kHz downsampled channels
            data = data.data;
        end
        
        DATA_MAT_EDF(ii_c+ii_cc,:) = data;
        EDFHeader.labels{ii_c+ii_cc} = sprintf('MACRO-%d-%s',channel,header.MontageCell{ii_c+ii_cc,2});
    end
    
    % add spike sorted cells to list
    if header.N_spikeSortedCells > 0
        Nchannels =  size(DATA_MAT_EDF,1);
        Nsamples = size(DATA_MAT_EDF,2);
        source_folder = header.processed_AverageRef;
        load(fullfile(source_folder,sprintf('%s_spike_timestamps_post_processing',header.id)))
        cnt = 0;
        for jj = [3 20 25 28 44 46]; % 1:length(micro_channels_spike_summary.spike_timestamps)
            cnt = cnt + 1;
            area_name = micro_channels_spike_summary.unit_list_xls(jj).Location;
            Cluster = micro_channels_spike_summary.unit_list_xls(jj).Cluster;
            Channel = micro_channels_spike_summary.unit_list_xls(jj).Channel;
            
            Unit_spike_timestamps_ms = micro_channels_spike_summary.spike_timestamps{jj};
            DATA_MAT_EDF(Nchannels+cnt,:) = zeros(Nsamples,1);
            DATA_MAT_EDF(Nchannels+cnt,int64(sort(unique([Unit_spike_timestamps_ms-1, Unit_spike_timestamps_ms,Unit_spike_timestamps_ms+1])))) = int8(100);
            EDFHeader.labels{Nchannels+cnt} = sprintf('%s-SU-Ch%dU%d',area_name,Channel,Cluster);
        end
        
    end
    
    %%% 10 Annotation
    % header.annotation.event     structure of cells whith event name
    % header.annotation.duration   vector with event duration (seconds)
    % header.annotation.starttime  vector with event startime  (seconds)
    Nchannels =  size(DATA_MAT_EDF,1);
    NStimCh = Nchannels + 1;
    Nsamples = size(DATA_MAT_EDF,2);
    A = zeros(Nsamples,1);
    EDFHeader.labels{NStimCh} = sprintf('stimulation');
    whatToDo2 = questdlg('Would you like to add stim channel?','EDF','YES','NO','NO');
    
    if strcmp(whatToDo2,'YES')
        
        whatToDo2 = questdlg('BR\NLX based stim?','EDF','BR','NLX','BR');
        if strcmp(whatToDo2,'BR')
            tstim = EXP_DATA.t_stim_TRAIN_END_MICRO_based_sec *1000; % msec, NLX time
            for jj_st = 1:length(tstim)
                A(int64(sort(unique([tstim-1, tstim,tstim+1])))) = 100;
            end
            DATA_MAT_EDF(NStimCh,:) = A;
            
        else
            
            if isfield(EXP_DATA,'stimTiming')
                tstim = EXP_DATA.stimTiming.validatedTTL_NLX;
            else
                TTL_struct = load(fullfile(header.processedDataPath,'STIM_TTLs_NLX.mat'),'TTL_struct'); TTL_struct = TTL_struct.TTL_struct;
                tstim = TTL_struct(1).t_stim_TRAIN_END_MICRO_based_sec *1000; % msec, NLX time
                tstim(isnan(tstim)) = [];
            end
            for jj_st = 1:length(tstim)
                A(int64(sort(unique([tstim-1, tstim,tstim+1])))) = 100;
            end
            DATA_MAT_EDF(NStimCh,:) = A;
        end
    end
    
    EDFHeader.samplerate = Consts.DOWNSAMPLED_FREQUENCY;
    EDFHeader.annotation = [];
    EDFHeader.annotation.event = [];
    EDFHeader.annotation.duration = [];
    EDFHeader.annotation.starttime = [];
    SaveEDF(fullfile(header.processed_MACRO,sprintf('%s_EXP%d_EDF_%s_1kHz.edf',header.id,header.experimentNum,str)), DATA_MAT_EDF', EDFHeader)
    
end

%%%--------------------------------------------------------------------------

whatToDo = questdlg('Extracting time from NLX files?','Q','YES','NO','NO');
if strcmp(whatToDo,'YES')
    DATA_DIR = inputdlg('Enter folder path:');
    DATA_DIR = DATA_DIR{1};
    CH = inputdlg('Enter channel name:');
    CH = CH{1};
    filelist = dir(fullfile(DATA_DIR,sprintf('%s*',CH)));
    FieldSelection = [0 0 0 0 1] ; % Will read all the variables from the file
    ExtractionMode = 2 ; % read only relevant session according to given timestamps
    ExtractionModeArray = [1 600000];
    
    for ii = 3:length(filelist)
        filename_CSC_EEG_in = fullfile(DATA_DIR, filelist(ii).name);
        a = dir(filename_CSC_EEG_in);
        if( a.bytes == 16384 )
            disp(sprintf('file %s is empty (%d bytes)',files(ii).name,files(ii).bytes))
            continue
        end
        
        [ ~, NlxHeader] = ...
            Nlx2MatCSC( filename_CSC_EEG_in, FieldSelection, 1, ExtractionMode, ExtractionModeArray ) ;
        
        filelist(ii).name
        NlxHeader(8)
        NlxHeader(9)
        keyboard
    end
    
end %function