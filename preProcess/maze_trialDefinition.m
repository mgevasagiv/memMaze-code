function pt_onsetPoints = maze_trialDefinition(ptname)

global subj_timestamps_folder

%% get event related 
% get decision points
timestamp_table = readtable(fullfile(subj_timestamps_folder, sprintf('%s_timestamp_table.xlsx',ptname)));

% ASK APRIL - this field changed into - Test_contextual_success_based_on_8_neighbors
count = 1;
for idx = 1:size(timestamp_table,1)
    if timestamp_table(idx,'Test_contextual_success').(1) == 1
        maze_set(count) = timestamp_table(idx,'Maze').(1);
        count = count + 1;
    end
end
good_maze_set = unique(maze_set); % These are the successful mazes in the contextual test

block_types = {'D','N','X','G'};
pt_onsetPoints.block_types = block_types;

%% for all mazes, each rep
for i = 1:length(block_types)
    block = block_types{i};
    for rep = 1:3
        
        clear all_D_idx
        count = 1;
        for idx = 1:size(timestamp_table,1)
            if (strcmp(cell2mat(timestamp_table(idx,'Type').(1)),block)) && (timestamp_table(idx,'Repetition').(1) == rep) % && (ismember(timestamp_table(idx,'Maze').(1), maze_set))
                all_D_idx(count) = idx;
                count = count + 1;
            end
        end
        
        clear onset_set 
        count = 1;
        for idx = 1:numel(all_D_idx)
            index = all_D_idx(idx);
            onset_i = timestamp_table(index,'Time_onset').(1);
            onset_set(count) = onset_i;
            count = count+1;
        end
        eval(['data_structure.',block,'_rep',num2str(rep),'= onset_set;']);
        
%         extracted_data = nan(numel(onset_set),totSamples,length(chan_nums));
%         for count = 1:numel(chan_nums)
%             for iTrial = 1:numel(onset_set)
%                 % extract data
%                 onset = onset_set(iTrial);
%                 onset_sample = (onset/1000000)*fsample;
%                 extracted_data(iTrial,:,count) = data_reshaped(onset_sample+plot_range(1):onset_sample+plot_range(2)-1,chan_nums(count));
%             end
%         end
%         eval([sub,'_data_structure.',block,'_rep',num2str(rep),'= extracted_data;']);
%         clear extracted_data all_D_idx  D_onset_set
    end
    
end

%% for good mazes, each rep
for i = 1:length(block_types)
    block = block_types{i};
    for rep = 1:3
        clear all_D_idx
        count = 1;
        for idx = 1:size(timestamp_table,1)
            if (strcmp(cell2mat(timestamp_table(idx,'Type').(1)),block)) && (timestamp_table(idx,'Repetition').(1) == rep)  && (ismember(timestamp_table(idx,'Maze').(1), good_maze_set))
                all_D_idx(count) = idx;
                count = count + 1;
            end
        end
        
        clear onset_set
        count = 1;
        for idx = 1:numel(all_D_idx)
            index = all_D_idx(idx);
            D_onset = timestamp_table(index,'Time_onset').(1);
            onset_set(count) = D_onset;
            count = count+1;
        end
                
        eval(['data_structure.',block,'_rep',num2str(rep),'_good_maze= onset_set;']);
        clear extracted_data all_D_idx  D_onset_set
    end
end

%% for rep2-3 & all mazes
for i = 1:length(block_types)
    block = block_types{i};
    
    
    count = 1;
    for idx = 1:size(timestamp_table,1)
        if (strcmp(cell2mat(timestamp_table(idx,'Type').(1)),block)) && ((timestamp_table(idx,'Repetition').(1) == 2) || (timestamp_table(idx,'Repetition').(1) == 3))
                all_D_idx(count) = idx;
                count = count + 1;
        end
    end
    
    count = 1;
    for idx = 1:numel(all_D_idx)
        index = all_D_idx(idx);
        D_onset = timestamp_table(index,'Time_onset').(1);
        onset_set(count) = D_onset;
        count = count+1;
    end
    
    eval(['data_structure.',block,'_rep2_3 = onset_set;']);
    clear  all_D_idx  onset_set
end

%% for rep2-3 & good mazes
for i = 1:length(block_types)
    block = block_types{i};
    
    
    count = 1;
    for idx = 1:size(timestamp_table,1)
        if (strcmp(cell2mat(timestamp_table(idx,'Type').(1)),block)) && ((timestamp_table(idx,'Repetition').(1) == 2) || (timestamp_table(idx,'Repetition').(1) == 3)) && (ismember(timestamp_table(idx,'Maze').(1), good_maze_set))
                all_D_idx(count) = idx;
                count = count + 1;
        end
    end
    
    count = 1;
    for idx = 1:numel(all_D_idx)
        index = all_D_idx(idx);
        D_onset = timestamp_table(index,'Time_onset').(1);
        onset_set(count) = D_onset;
        count = count+1;
    end
        
    eval(['data_structure.',block,'_rep2_3_good_maze = onset_set;']);
    clear  all_D_idx  onset_set

end

pt_onsetPoints.onsetPoints = data_structure;
