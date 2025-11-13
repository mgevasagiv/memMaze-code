function generate_data_structure_goalDecode(sub,subj_folder,chan_nums,plot_range,exp_num)
% sub: subject number
% subj_folder: input data folder for the subject
% chan_nums: channel IDs intended to obtain for the data structure
% plot_range: time sample range
MACRO_MONTAGE_FOLDER = 'E:\MAZE\data_p\MACRO_MONTAGE\';
BEHAV_TRIALS_FOLDER = 'E:\MAZE\data_p\populationAnalysis\MACRO\BEHAV_TRIALS';
addpath( 'C:\Users\mgeva\Documents\GitHub\external\fieldtrip-20230613\fieldtrip-20230613\' )
addpath( 'C:\Users\mgeva\Documents\GitHub\external\fieldtrip-20230613\fieldtrip-20230613\fileio\' )
addpath(genpath('C:\Users\mgeva\Documents\GitHub\closedLoop-code\tools\shadedErrorbar\'))
addpath('C:\Users\mgeva\Documents\GitHub')
addpath(genpath('C:\Users\mgeva\Documents\GitHub\memMaze-code'))
addpath(genpath('C:\Users\mgeva\Documents\GitHub\closedLoop-code\epilepticActivity_IEEG-code'))
rmFieldTripAccessPath

%% load all channel data
load(fullfile(MACRO_MONTAGE_FOLDER, sprintf('p%s',sub),sprintf('EXP%s',exp_num), 'MacroMontage.mat'));
numChan = length(MacroMontage);
chan_nums  = getDeepMacros(MacroMontage); % get deep contact on each electrode

rmv_id = [];
for n = 1:length(chan_nums)
    if isempty(dir(fullfile(subj_folder, sprintf('CSC%d.mat',chan_nums(n)))))
        rmv_id = [rmv_id ,n];
    end
end
chan_nums(rmv_id) = [];

if isempty(subj_folder)
    subj_folder = header.processed_MACRO;
end
file = load(fullfile(subj_folder, 'CSC1.mat'));
fsample = file.CSC_Sampling_Rate_Hz;
totSamples = diff(plot_range);
if fsample ~= 1000; error('totSamples should be adjusted for different sampling rates'); end
IIS_det = SpikeWaveDetector;
IIS_det.samplingRate = fsample;
IIS_det.useAmp = 1;
IIS_det.SDthresholdAmp = 5;
IIS_det.SDthresholdEnv = 3;
IIS_det.blockSizeSec = 120;
IIS_det.useConjAmpEnv = 1;

data = file.data;
data= data(:);
sEEG_data = nan(size(data,1),length(chan_nums));
% PreProcessing to generate trials:
% -- 60Hz was removed on iEEG extraction [see remove_line_noise()]
% -- bandpass to remove high-freq noise [1 250]
% -- interictal spike detection, followed by removal of IED events
% 
for n = 1:length(chan_nums)
    try
        file = load(fullfile(subj_folder, sprintf('CSC%d.mat',chan_nums(n))));
    catch
        fprintf('CSC%d.mat is missing\n',chan_nums(n))
        continue
    end
    data = file.data;
    low_cut = 1;
    high_cut = 250;
    data_lp = bandpass(data, fsample, low_cut, high_cut); % Low pass siganl
    [peakTimes, peakStats] = IIS_det.detectTimes(data(:)', true);
    disp(sprintf('removing %d spikes from ch %d',length(peakTimes),chan_nums(n)))
    data_lp_clean = data_lp;
    data_lp_clean(peakTimes-50: peakTimes+50) = 0; % remove detected interictal spikes 

    if length(data_lp_clean) ~= length(data)
        error('CSC - wrong length')
    else
        sEEG_data(:,n) = data_lp_clean(:)';
    end
end

% get unique electrode shank name
prev_str = MacroMontage(1).Area;
electrode_names{1} = MacroMontage(1).Area;
num_shank = 2;
num_electrodes = 0;
for i = 1:length(MacroMontage)
    electrode_name = MacroMontage(i).Area;
    if ~strcmp(electrode_name,prev_str)
        prev_str = electrode_name;
        electrode_names{num_shank} = electrode_name;
        electrode_num{num_shank-1} = num_electrodes;
        num_shank = num_shank + 1;
        num_electrodes = 0;
    end
    num_electrodes = num_electrodes + 1;
end
electrode_num{num_shank-1} = num_electrodes;


%% load channel info
count = 1;
for i = 1:length(electrode_names)
    for j = 1:cell2mat(electrode_num(i))
        name = [cell2mat(electrode_names(i)),num2str(j)];
        channel_data{count} = name;
        count = count + 1;
    end
end
channel_data = channel_data(chan_nums)';

for i = 1:length(channel_data)
    if strfind(channel_data{i}, 'WM')
        rmv_id = [rmv_id i];
    end
end
% RMV white matter electrodes from analysis
chan_nums(rmv_id) = [];
channel_data(rmv_id)= [];

% preprocess all data
cfg = [];
cfg.lpfilter      = 'yes';
cfg.lpfreq        = 200;
combined_data.label = channel_data;
combined_data.trial = sEEG_data';
combined_data.fsample = file.CSC_Sampling_Rate_Hz;
combined_data.time = [0 : size(data,1)-1 ]./ combined_data.fsample;
% combined_data               = ft_preprocessing(cfg,combined_data);

%% get event related
% get decision points
timestamp_folder = 'E:\MAZE\ptData';
if str2num(exp_num) == 2
    timestamp_table = readtable(fullfile(timestamp_folder, [sub '_e2_timestamp_table.csv']));
    
else
    try
        timestamp_table = readtable(fullfile(timestamp_folder, [sub '_timestamp_table.csv']));
    catch
        timestamp_table = readtable(fullfile(timestamp_folder, [sub '_timestamp_table.xlsx']));
    end
end

try
    mm = matfile(fullfile(BEHAV_TRIALS_FOLDER, [sub,'_',exp_num,'behav_data_structure_goal_decode.mat']));
    behav_data_structure = mm.behav_data_structure;
    mazeCompletionTable = mm.mazeCompletionTable;
catch
    error('maze info is missing')
end
rmv_idx = find(timestamp_table.Time_onset <0 );
timestamp_table.Repetition(rmv_idx) = nan;

%% finding good mazes based on surrounding 8 mazes
count = 1;
for idx = 1:size(timestamp_table,1)
    if timestamp_table(idx,'Test_contextual_success').(1) == 1 % Switch to 'Test_contextual_success_based_on_8_neighbors'
        maze_set(count) = timestamp_table(idx,'Maze').(1);
        count = count + 1;
    end
end
good_maze_set = unique(maze_set);

% Focus on exploration trials
block_types = {'X','N','G'};
uSec2Sec = 1e6;
%% for all mazes, each rep
explore_data_structure.pt_ID = sub;


for i = 1:length(block_types)
    block = block_types{i};
    clear all_D_idx state_idx_good_maze D_onset_set D_maze_id
     
    for rep = 1
        
        % First collect all relevant events
        count = 1;
        for idx = 1:size(timestamp_table,1)
            if (strcmp(cell2mat(timestamp_table(idx,'Type').(1)),block)) && (timestamp_table(idx,'Repetition').(1) == rep) % && (ismember(timestamp_table(idx,'Maze').(1), maze_set))
                if (timestamp_table(idx,'Time_onset').(1)/uSec2Sec)*fsample+plot_range(1) > 0
                    all_D_idx(count) = idx;
                    count = count + 1;
                end
            end
        end
        
        % second -  focus on those that are part of a succesfully remembered goal-maze 
        count = 1;
        for idx = 1:size(timestamp_table,1)
            if (strcmp(cell2mat(timestamp_table(idx,'Type').(1)),block)) && (timestamp_table(idx,'Repetition').(1) == rep)  && (ismember(timestamp_table(idx,'Maze').(1), good_maze_set))
                if (timestamp_table(idx,'Time_onset').(1)/uSec2Sec)*fsample+plot_range(1) > 0
                    state_idx_good_maze(count) = idx;
                    count = count + 1;
                end
            end
        end
        trials_good_maze = find(ismember(all_D_idx, state_idx_good_maze)); % indices of good mazes 

        count = 1;
        for idx = 1:numel(all_D_idx)
            index = all_D_idx(idx);
            D_onset_set(count) = timestamp_table(index,'Time_onset').(1);
            D_maze_id(count) = timestamp_table(index,'Maze').(1);
            count = count+1;
        end
        
        onset_set = D_onset_set; % Full set of behavioral points
        q_loc = mazeCompletionTable.goal_quadrant(D_maze_id + 1); % mazes numbering starts at 0;

        if strcmp(block,'N')
            neutral_set = floor((onset_set/uSec2Sec)*fsample);
        end
        
        extracted_data = nan(numel(onset_set),totSamples,length(chan_nums));
        zsore_threshold = 7;
        for count = 1:numel(chan_nums)
            zscore_data = zscore(combined_data.trial(count,:));
            for iTrial = 1:numel(onset_set)
                cnt = 1;
                % extract data
                onset = onset_set(iTrial);
                onset_sample = floor((onset/uSec2Sec)*fsample);
                trial_indices = onset_sample+plot_range(1):onset_sample+plot_range(2)-1;
                if sum(abs(zscore_data(trial_indices)) > zsore_threshold)
                    cnt = cnt + 1;                   
                    extracted_data(iTrial,:,count) = zeros(1,length(trial_indices));
                else
                    extracted_data(iTrial,:,count) = combined_data.trial(count,trial_indices);
                end
            end
            disp(sprintf(' %d trial removed',cnt))
        end

        eval(['explore_data_structure.',block,'_rep',num2str(rep),'= extracted_data;']);
        eval(['explore_data_structure.',block,'_rep',num2str(rep),'_good_maze_trials','= trials_good_maze;']);
        eval(['explore_data_structure.',block,'_rep',num2str(rep),'_maze_id_set','= D_maze_id;']);   
        eval(['explore_data_structure.',block,'_rep',num2str(rep),'_maze_trials_goal_q','= q_loc;']);

        % sanity check

        if (length(D_maze_id) ~= length(q_loc))
            error('vector size mismatch')
        end

        clear extracted_data all_D_idx  D_onset_set state_idx_good_maze trials_good_maze


    end
end

experienced_data_structure.pt_ID = sub;

for i = 1:length(block_types)

    block = block_types{i};
    clear all_D_idx state_idx_good_maze D_onset_set D_maze_id

    for rep = 2:3
        
        % First collect all relevant events
        count = 1;
        for idx = 1:size(timestamp_table,1)
            if (strcmp(cell2mat(timestamp_table(idx,'Type').(1)),block)) && (timestamp_table(idx,'Repetition').(1) == rep) % && (ismember(timestamp_table(idx,'Maze').(1), maze_set))
                if (timestamp_table(idx,'Time_onset').(1)/uSec2Sec)*fsample+plot_range(1) > 0
                    all_D_idx(count) = idx;
                    count = count + 1;
                end
            end
        end
        
        % second -  focus on those that are part of a succesfully remembered goal-maze 
        count = 1;
        for idx = 1:size(timestamp_table,1)
            if (strcmp(cell2mat(timestamp_table(idx,'Type').(1)),block)) && (timestamp_table(idx,'Repetition').(1) == rep)  && (ismember(timestamp_table(idx,'Maze').(1), good_maze_set))
                if (timestamp_table(idx,'Time_onset').(1)/uSec2Sec)*fsample+plot_range(1) > 0
                    state_idx_good_maze(count) = idx;
                    count = count + 1;
                end
            end
        end
        trials_good_maze = find(ismember(all_D_idx, state_idx_good_maze)); % indices of good mazes 
        
        count = 1;
        for idx = 1:numel(all_D_idx)
            index = all_D_idx(idx);
            D_onset = timestamp_table(index,'Time_onset').(1);
            D_onset_set(count) = D_onset;
            D_maze_id(count) = timestamp_table(index,'Maze').(1);
            count = count+1;
        end
        
        onset_set = D_onset_set; % Full set of behavioral points
                q_loc = mazeCompletionTable.goal_quadrant(D_maze_id + 1); % mazes numbering starts at 0;

        if strcmp(block,'N')
            neutral_set = floor((onset_set/uSec2Sec)*fsample);
        end
        
        extracted_data = nan(numel(onset_set),totSamples,length(chan_nums));
        zsore_threshold = 7;
        for count = 1:numel(chan_nums)
            zscore_data = zscore(combined_data.trial(count,:));
            for iTrial = 1:numel(onset_set)
                cnt = 1;
                % extract data
                onset = onset_set(iTrial);
                onset_sample = floor((onset/uSec2Sec)*fsample);
                trial_indices = onset_sample+plot_range(1):onset_sample+plot_range(2)-1;
                if sum(abs(zscore_data(trial_indices)) > zsore_threshold)
                    cnt = cnt + 1;                   
                    extracted_data(iTrial,:,count) = zeros(1,length(trial_indices));
                else
                    extracted_data(iTrial,:,count) = combined_data.trial(count,trial_indices);
                end
            end
            disp(sprintf(' %d trial removed',cnt))
        end

        eval(['experienced_data_structure.',block,'_rep',num2str(rep),'= extracted_data;']);
        eval(['experienced_data_structure.',block,'_rep',num2str(rep),'_good_maze_trials','= trials_good_maze;']);
        eval(['experienced_data_structure.',block,'_rep',num2str(rep),'_maze_id_set','= D_maze_id;']);
        eval(['experienced_data_structure.',block,'_rep',num2str(rep),'_maze_trials_goal_q','= q_loc;']);


        % sanity check
        if (length(D_maze_id) ~= length(q_loc))
            error('vector size mismatch') 
        end

        clear extracted_data all_D_idx  D_onset_set state_idx_good_maze trials_good_maze D_maze_id

    end
end


eval([sub,'_data_structure.pt_ID = sub;']);
eval([sub,'_data_structure.sampling_rate = fsample;']);
eval([sub,'_data_structure.channels_included = chan_nums;']);
eval([sub,'_data_structure.channels_included_areas = channel_data;']);
eval([sub,'_data_structure.experiment_num =', num2str(exp_num),';']);

save(fullfile(BEHAV_TRIALS_FOLDER, [sub,'_',exp_num,'_data_structure_goalDecode.mat']),[sub,'_data_structure'],'explore_data_structure','experienced_data_structure');

end


%% help functions

function BP = bandpass(timecourse, SamplingRate, low_cut, high_cut, filterOrder)
%constants for bandpass
defaultFilterOrder = 1;
nanWarning = 0.01;

%bandpass code - from Maya
if (nargin < 6)
    filterOrder = defaultFilterOrder;
end

% Maya GS - handle NAN values
indices = find(isnan(timecourse));
if length(indices) > nanWarning*length(timecourse)
    warning('many NaN values in filtered signal')
end
timecourse(indices) = 0;
%
if high_cut == inf
    [b, a] = butter(filterOrder,(low_cut/SamplingRate)*2,'high');
else
    [b, a] = butter(filterOrder, [(low_cut/SamplingRate)*2 (high_cut/SamplingRate)*2]);
end
BP = filtfilt(b, a, timecourse );
BP(indices) = NaN;
end