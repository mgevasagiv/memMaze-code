function generateBehavDataStructure(sub,subj_folder,exp_num)
% sub: subject number
% subj_folder: input data folder for the subject
% chan_nums: channel IDs intended to obtain for the data structure
% plot_range: time sample range
MACRO_MONTAGE_FOLDER = 'E:\MAZE\data_p\MACRO_MONTAGE\';
BEHAV_TRIALS_FOLDER = 'E:\MAZE\data_p\populationAnalysis\MACRO\BEHAV_TRIALS';
addpath( 'C:\Users\mgeva\Documents\GitHub\external\fieldtrip-20230613\fieldtrip-20230613\' )
addpath( 'C:\Users\mgeva\Documents\GitHub\external\fieldtrip-20230613\fieldtrip-20230613\fileio\' )
addpath(genpath('C:\Users\mgeva\Documents\GitHub\closedLoop-code\tools\shadedErrorbar\'))
addpath('C:\Users\mgeva\Documents\GitHub\')
rmFieldTripAccessPath

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

rmv_idx = find(timestamp_table.Time_onset <0 );
timestamp_table.Repetition(rmv_idx) = nan;

%%
% Look at reaction times
block_types = {'D','X','G','N',};
uSec2Sec = 1e6;
fsample = 1000;
plot_range= 500;
%% for all mazes, each rep
for i = 1:length(block_types)
    block = block_types{i};
    for rep = 1:3
        maze_id_found = nan(1,24);
        % First collect all relevant events
        count = 1;
        for idx = 1:size(timestamp_table,1)
            if (strcmp(cell2mat(timestamp_table(idx,'Type').(1)),block)) && (timestamp_table(idx,'Repetition').(1) == rep) && (timestamp_table.Block(idx) < 10) % looking at first 3 blocks
                if (timestamp_table(idx,'Time_onset').(1)/uSec2Sec)*fsample+plot_range(1) > 0
                    if ismember(i,2:3)
                        if isnan(timestamp_table.Maze(idx))
                            continue
                        else
                            maze_id_found(timestamp_table.Maze(idx)+1) = 1;
                        end
                    end
                    all_idx(count) = idx;
                    reaction_times(count) = timestamp_table.Time_offset(idx)-timestamp_table.Time_onset(idx);
                    maze_id(count) = timestamp_table.Maze(idx);
                    count = count + 1;
                end
            end
        end
        eval(['behav_data_structure.',block,'_rep',num2str(rep),'= all_idx;']);
        eval(['behav_data_structure.',block,'_rep',num2str(rep),'_reaction_times','= reaction_times;']);
        eval(['behav_data_structure.',block,'_rep',num2str(rep),'_maze_id','= maze_id;']);
        clear all_idx reaction_times maze_id
    end
    
end

if (sum(~ismember([0:23],behav_data_structure.G_rep3_maze_id)) || sum(~ismember([0:23],behav_data_structure.G_rep2_maze_id)) || sum(~ismember([0:23],behav_data_structure.G_rep1_maze_id)))
    warning('maze is missing')
end

% Prep a table with time to complete each maze

% Look at path length to goal
%% for all mazes, each rep
mazeIds = [0:23]';
for rep = 1:3
    for maze_id = 1:length(mazeIds)
        % First collect all relevant events
        % search for start and end times for each maze/rep
        block_types = {'X','G'};
        count = 1;
        all_idx = [];
        for i = 1:length(block_types)
            block = block_types{i};
            for idx = 1:size(timestamp_table,1)
                if (strcmp(cell2mat(timestamp_table(idx,'Type').(1)),block)) && (timestamp_table(idx,'Repetition').(1) == rep) && (timestamp_table(idx,'Maze').(1)== mazeIds(maze_id))
                    if (timestamp_table(idx,'Time_onset').(1)/uSec2Sec)*fsample+plot_range(1) > 0
                        all_idx(count) = idx;
                        maze_id_validation(maze_id,rep) = timestamp_table.Maze(idx);
                        count = count + 1;
                    end
                end
            end
        end
        if isempty(all_idx) || (length(all_idx) < 2)
            mazeCompletion(maze_id,rep) = NaN;
        else
            mazeCompletion(maze_id,rep) = timestamp_table.Time_offset(all_idx(2))-timestamp_table.Time_onset(all_idx(1));
        end
    end
end

mazeCompletionTable = table( mazeIds, mazeCompletion(:,1),mazeCompletion(:,2),mazeCompletion(:,3),'VariableNames',{'maze_id','solve_ms_rep1','solve_ms_rep2','solve_ms_rep3'});

save(fullfile(BEHAV_TRIALS_FOLDER, [sub,'_',exp_num,'behav_data_structure.mat']),'behav_data_structure','mazeCompletionTable');

end