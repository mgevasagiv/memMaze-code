function generateBehavDataStructureGoalDecoding(sub,subj_folder,exp_num)
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

try
    goal_loc_file = 'E:\MAZE\ptData\maze_goals_with_quadrants.csv';
    goal_maze_table = readtable(goal_loc_file);
catch
    error('maze-goal pairing not found')
end

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
        block_types = {'X','E'};
        count = 1;
        all_idx = [];
        for i = 1:length(block_types)
            block = block_types{i};
            for idx = 1:size(timestamp_table,1)
                if ( strcmp(cell2mat(timestamp_table(idx,'Type').(1)),block) && (timestamp_table(idx,'Repetition').(1) == rep) && ...
                        (timestamp_table(idx,'Maze').(1)== mazeIds(maze_id)) && (timestamp_table(idx,'Block').(1) <= 9) )
                    if (timestamp_table(idx,'Time_onset').(1)/uSec2Sec)*fsample+plot_range(1) > 0
                        all_idx(count) = idx(1); % sometimes people go back accidently before the maze switches 
                        maze_id_validation(maze_id,rep) = timestamp_table.Maze(idx(1));
                        count = count + 1;
                    end
                end
            end
        end
        if isempty(all_idx) || (length(all_idx) < 2)
            mazeCompletion(maze_id,rep) = NaN;
            steps_to_completion(maze_id,rep) = NaN;

            continue
        else
            mazeCompletion(maze_id,rep) = timestamp_table.Time_offset(all_idx(2))-timestamp_table.Time_onset(all_idx(1));
        end
        
        
        if ~strcmp(timestamp_table(all_idx(1),'Type').(1),'X') || ~strcmp(timestamp_table(all_idx(2),'Type').(1),'E')
            error('bad trial segmentation'); end

        steps_to_completion(maze_id,rep) = all_idx(2)-all_idx(1);

    end
end
maze_goal_Q = goal_maze_table.quadrant(1:24);

clear extracted_data

str = 'test';
testing_indices = find(~isnan(timestamp_table.Test_contextual_success) & timestamp_table.Block <= 9 );
corr_indices = testing_indices(find(timestamp_table.Test_contextual_success(testing_indices) == 1));

Repitition_correct  =  timestamp_table.Repetition(corr_indices);
Maze_correct =  timestamp_table.Maze(corr_indices);
inCorrectTest = timestamp_table.Test_contextual_success(testing_indices) == 0;


Mat_Maze_sol = zeros(24,3);
for ii = 1:length(Maze_correct)
    Mat_Maze_sol(Maze_correct(ii)+1,Repitition_correct(ii)) = 1;
end

mazeCompletionTable = table( mazeIds, maze_goal_Q, mazeCompletion(:,1),mazeCompletion(:,2),mazeCompletion(:,3),...
    steps_to_completion(:,1),steps_to_completion(:,2),steps_to_completion(:,3),...
    'VariableNames',{'maze_id','goal_quadrant','solve_ms_rep1','solve_ms_rep2','solve_ms_rep3',...
        'steps_to_solve_ms_rep1','steps_to_solve_ms_rep2','steps_to_solve_ms_rep3'});

% prep figure
newA4figure([sub,'_',exp_num,'_Trialtime'])
axes('Position',[0.1 0.1 0.3 0.4])
hold all
mazeCompletion_s = mazeCompletion(:,1:3)./mazeCompletion(:,1);



plot(mazeCompletion_s','.--','color',[0.5 0.5 0.5])
my_boxplot(mazeCompletion_s(:,1),1,'k',[12 0.5])
my_boxplot(mazeCompletion_s(:,2),2,'k',[12 0.5])
my_boxplot(mazeCompletion_s(:,3),3,'k',[12 0.5])
set(gca,'xtick',[1:3])
xlabel('Repitition')
ylabel('Trial time (norm)')
saveas(gcf, fullfile('E:\MAZE\data_p\populationAnalysis\behavFigures',[get(gcf,'name'),'.tif']))

axes('Position',[0.5 0.1 0.3 0.4])
hold all
% steps_to_completion = steps_to_completion(:,1:3)./steps_to_completion(:,1);

plot(steps_to_completion','.--','color',[0.5 0.5 0.5])
my_boxplot(steps_to_completion(:,1),1,'k',[12 0.5])
my_boxplot(steps_to_completion(:,2),2,'k',[12 0.5])
my_boxplot(steps_to_completion(:,3),3,'k',[12 0.5])

for ii = 1:length(Mat_Maze_sol)
    for jj = 1:3
        if(Mat_Maze_sol(ii,jj) == 1);
            plot(jj,steps_to_completion(ii,jj),'go','MarkerSize',10,'LineWidth',2.5)
        else
            plot(jj,steps_to_completion(ii,jj),'rd','MarkerSize',10,'LineWidth',2.5)
        end
    end
end

set(gca,'xtick',[1:3])
xlabel('Repitition')
ylabel('Steps per trial')

saveas(gcf, fullfile('E:\MAZE\data_p\populationAnalysis\behavFigures',[get(gcf,'name'),'.tif']))


save(fullfile(BEHAV_TRIALS_FOLDER, [sub,'_',exp_num,'behav_data_structure_goal_decode.mat']),'behav_data_structure','mazeCompletionTable');

save(fullfile(BEHAV_TRIALS_FOLDER, [sub,'_',exp_num,'behav_data_structure_goal_decode.mat']),'behav_data_structure','mazeCompletionTable');

end