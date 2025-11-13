% Analyze MACRO LFP during MAZE goal-oriented activity
clear all; close all

codePath_MAZE();
global patients expNames spikeSortingFolderAll data_p_path BEHAV_TRIALS_FOLDER MACROLFP_spectralAnalysisFolder

HPC_channels_spectralAnalysis = {[],... %17
    [],... % 18
    [],... % 19_1
    [1,103],... % 19_2 - Am, HPC
    [],... % 22 - no MTL
    [1,2,17,18],...%  23 - Am, HPC
    [27:29.37:39],... % 26, Hip, AM
    [1,2,13,14,23,24,35,36,47,48, 59,60],... % ir103, LAM, RAM, LHeadH, RHeadH, LTailH, RTailH
    [17 18 ],... % p30
    [1,2],... % p35
    [17:18,27:28.37:38,132:133,142:143,152:153],... % ir106
    }; 

OFC_channels_spectralAnalysis = {[],... % 17
    [],... % 18
    [],... % 19_1
    [57,80],... % 19_2
    [],... % 22
    [119,120],...% 23
    [47:49],... % 26
    [215, 227],... % ir103 left anterior cingulate, right anterior cingulate
    [68 69 122 123 126],... % p30
    [56,57,51,61],... % p35 OF and ACC
    [47:48,57:58,162:163,171:172],... % ir106
    };  

data_p_path = 'E:\MAZE\Data_p\';
spikeSortingFolderAll = fullfile(data_p_path,'singleUnit');

runData = [];
nPatients = length(patients);
for iPatient = 1:nPatients % Until Data from  I106 arrives
    runData(iPatient).patientName = patients{iPatient};
    runData(iPatient).EXP = num2str(expNames{iPatient}(4:end));
    
    runData(iPatient).DataFolder = [data_p_path, patients{iPatient},'\',expNames{iPatient},'\Denoised_Downsampled_InMicroVolt\MACRO'];
    
    %The folder+filename into which the spikes results is going to be stored or is
    %already stored if the spikes detection was already run (the folder should
    %exist)
    runData(iPatient).SpikesFileNames = fullfile(runData(iPatient).DataFolder, 'spikesResults',...
        sprintf('SpikesResults'));
    
    runData(iPatient).ExpDataFileName = [data_p_path,patients{iPatient},'\',expNames{iPatient},'\',patients{iPatient},'_',expNames{iPatient},'_dataset.mat'];
    
    runData(iPatient).microMontageFileName = [data_p_path,'MACRO_MONTAGE','\',patients{iPatient},'\',expNames{iPatient},'\microMontage.mat'];
    runData(iPatient).macroMontageFileName = [data_p_path,'MACRO_MONTAGE','\',patients{iPatient},'\',expNames{iPatient},'\MacroMontage.mat'];
    runData(iPatient).spikeData = [data_p_path,patients{iPatient},'\',expNames{iPatient},'\spikeSorting\',patients{iPatient},'_spike_timestamps_post_processing.mat'];
    
    runData(iPatient).OFC_channels_spectralAnalysis = OFC_channels_spectralAnalysis{iPatient};
    runData(iPatient).HPC_channels_spectralAnalysis = HPC_channels_spectralAnalysis{iPatient};

    % define 'legal' channels for analysis
    % remove pathological channnels
    % remove white matter channels
    runData(iPatient).channelsForAnalysis = [];
    results_folder_top = 'E:\Data_p\ClosedLoopDataset\';

    %if isempty(dir(fullfile(BEHAV_TRIALS_FOLDER,sprintf('*%s_%d*behav_data_structure_goal_decode.mat',runData(iPatient).patientName(4:end),str2num(expNames{iPatient}(end))))))
    %    generateBehavDataStructureGoalDecoding(runData(iPatient).patientName(2:end),runData(iPatient).DataFolder,runData(iPatient).EXP)
        plot_range_ms = [-500 2000];
        generate_data_structure_goalDecode(runData(iPatient).patientName(2:end) ,runData(iPatient).DataFolder,'',plot_range_ms, runData(iPatient).EXP )
    %end

    if isempty(dir(fullfile(BEHAV_TRIALS_FOLDER,sprintf('*%s_%d*data_structure.mat',runData(iPatient).patientName(4:end),str2num(expNames{iPatient}(end)))))) 
        generateBehavDataStructure(runData(iPatient).patientName(2:end),runData(iPatient).DataFolder,runData(iPatient).EXP)
    end
    
    if isempty(dir(fullfile(BEHAV_TRIALS_FOLDER,sprintf('*%s_%d*data_structure.mat',runData(iPatient).patientName(4:end),str2num(expNames{iPatient}(end))))))
        disp(sprintf('creating data structure for %s',patients{iPatient}))
        plot_range_ms = [-500 500];
        generate_data_structure(runData(iPatient).patientName(2:end) ,runData(iPatient).DataFolder,'',plot_range_ms, runData(iPatient).EXP )
    end

end

AMNA = AnalyzeMazeNeuralActivity;

% Look at MACRO activity
outputFileFolder = fullfile(data_p_path,'populationAnalysis\macroLFP_ERP');
outputFigureFolder = fullfile(data_p_path,'populationAnalysis\macroLFP_ERP\figures');

disp('G')
whatToRun.timeStampsType = 'G';
AMNA.MACROLFP_collectPerMazeNeuralMarkers( runData, outputFileFolder, outputFigureFolder, whatToRun);
disp('X')
whatToRun.timeStampsType = 'X';
AMNA.MACROLFP_collectPerMazeNeuralMarkers( runData, outputFileFolder, outputFigureFolder, whatToRun);
disp('D')
whatToRun.timeStampsType = 'D';
AMNA.MACROLFP_collectPerMazeNeuralMarkers( runData, outputFileFolder, outputFigureFolder, whatToRun);

% analyse the data population 
genDataFilesModeling() % aggregate files generate from 'MACROLFP_collectPerMazeNeuralMarkers'
useNeuralBehavTables()
genFigureGoalPlanningChanges()


% Spectral analysis 
% Uses mtspecgramc to calculate power spectrum for trials and for control
% trials 

outputFileFolder = fullfile(data_p_path,'populationAnalysis\macroLFP_spectralAnalysis');
outputFigureFolder = fullfile(data_p_path,'populationAnalysis\macroLFP_spectralAnalysis\figures');

whatToRun.timeStampsType = 'D';
whatToRun.singlePlots = false;
whatToRun.newtimef = true;
AMNA.MACROLFP_spectralAnalysis( runData, outputFileFolder, outputFigureFolder, whatToRun);

whatToRun.timeStampsType = 'X';
AMNA.MACROLFP_spectralAnalysis( runData, outputFileFolder, outputFigureFolder, whatToRun);

whatToRun.timeStampsType = 'G';
AMNA.MACROLFP_spectralAnalysis( runData, outputFileFolder, outputFigureFolder, whatToRun);

whatToRun.timeStampsType = 'Gm';
AMNA.MACROLFP_spectralAnalysis( runData, outputFileFolder, outputFigureFolder, whatToRun);

whatToRun.timeStampsType = 'Xm';
AMNA.MACROLFP_spectralAnalysis( runData, outputFileFolder, outputFigureFolder, whatToRun);

% whatToRun.timeStampsType = 'Test'; % NEED TO LOCATE TEST TIMESTAMPS
% AMNA.MACROLFP_spectralAnalysis( runData, outputFileFolder, outputFigureFolder, whatToRun);


outputFileFolder = fullfile(data_p_path,'populationAnalysis\macroLFP_spectralAnalysis');
AMNA.AMNA_collect_population_newtimeF(runData, outputFileFolder, outputFigureFolder);

outputFigureFolder = fullfile(data_p_path,'populationAnalysis\macroLFP_spectralAnalysis\figures\ERSP');
AMNA.ersp_stats(outputFileFolder, outputFigureFolder)

whatToRun.timeStampsType = 'G';
whatToRun.singlePlots = true;
whatToRun.isMTLPt = 1;
AMNA.MACROLFP_crossChannel_spectralAnalysis(runData, outputFileFolder, outputFigureFolder, whatToRun)


% Combine all data files
outputFileFolder = fullfile(data_p_path,'populationAnalysis\macroLFP_spectralAnalysis');
outputFigureFolder = fullfile(data_p_path,'populationAnalysis\macroLFP_spectralAnalysis\figures');

whatToRun.timeStampsType = 'D';
AMNA.AMNA_population_analysis(runData, outputFileFolder, outputFigureFolder, whatToRun);
whatToRun.timeStampsType = 'X';
AMNA.AMNA_population_analysis(runData, outputFileFolder, outputFigureFolder, whatToRun);
whatToRun.timeStampsType = 'G';
AMNA.AMNA_population_analysis(runData, outputFileFolder, outputFigureFolder, whatToRun);
whatToRun.timeStampsType = 'Gm';
AMNA.AMNA_population_analysis(runData, outputFileFolder, outputFigureFolder, whatToRun);
whatToRun.timeStampsType = 'Xm';
AMNA.AMNA_population_analysis(runData, outputFileFolder, outputFigureFolder, whatToRun);

% Single plots per subject - ERPs and stats
whatToRun.singlePlots = true;
whatToRun.plotSig = false;
whatToRun.timeStampsType = 'G';
AMNA.AMNA_population_figures(runData, outputFileFolder, outputFigureFolder, whatToRun)
whatToRun.timeStampsType = 'Gm';
AMNA.AMNA_population_figures(runData, outputFileFolder, outputFigureFolder, whatToRun)
whatToRun.timeStampsType = 'Xm';
AMNA.AMNA_population_figures(runData, outputFileFolder, outputFigureFolder, whatToRun)
whatToRun.timeStampsType = 'D';
AMNA.AMNA_population_figures(runData, outputFileFolder, outputFigureFolder, whatToRun)

plotCoherencePop.m