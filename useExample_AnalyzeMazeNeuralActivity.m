global patients expNames spikeSortingFolderAll data_p_path 
patients = {'pda017','pda019'};
expNames = {'EXP1','EXP1'};

microChannels_spectralAnalysis = {[2:11],...
                                   [1:16]}; 

data_p_path = 'E:\MAZE\Data_p\';
spikeSortingFolderAll = fullfile(data_p_path,'singleUnit');

runData = [];
nPatients = length(patients);
for iPatient = 1:nPatients
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
    
    runData(iPatient).microChannels_spectralAnalysis = microChannels_spectralAnalysis{iPatient};
    
    results_folder_top = 'E:\Data_p\ClosedLoopDataset\';

end

AMNA = AnalyzeMazeNeuralActivity;
 
%% 

% Look at single unit activity - ERP
for iiP = 1:length(runData)
    
    figureFolder = fullfile(data_p_path,'populationAnalysis\neuralUnits_ERP_figures');
    whatToRun.timeStampsType = 'D';
    fileNameResults = fullfile(data_p_path, 'populationAnalysis\neuralUnits_trialData', sprintf('neuralUnitsTrials_%s_EXP%s_cond%s',runData(iiP).patientName,runData(iiP).EXP,whatToRun.timeStampsType) );
    results = AMNA.getStimulusTriggeredFireRateAllAreas(runData(iiP), fileNameResults, whatToRun);
    PrintActiveFigs(figureFolder)
    
    whatToRun.timeStampsType = 'X';
    fileNameResults = fullfile(data_p_path, 'populationAnalysis\neuralUnits_trialData', sprintf('neuralUnitsTrials_%s_EXP%s_cond%s',runData(iiP).patientName,runData(iiP).EXP,whatToRun.timeStampsType) );
    results = AMNA.getStimulusTriggeredFireRateAllAreas(runData(iiP), fileNameResults, whatToRun);
    PrintActiveFigs(figureFolder)
    
    whatToRun.timeStampsType = 'G';
    fileNameResults = fullfile(data_p_path, 'populationAnalysis\neuralUnits_trialData', sprintf('neuralUnitsTrials_%s_EXP%s_cond%s',runData(iiP).patientName,runData(iiP).EXP,whatToRun.timeStampsType) );
    results = AMNA.getStimulusTriggeredFireRateAllAreas(runData(iiP), fileNameResults, whatToRun);
    PrintActiveFigs(figureFolder)
    
    whatToRun.timeStampsType = 'Gm';
    fileNameResults = fullfile(data_p_path, 'populationAnalysis\neuralUnits_trialData', sprintf('neuralUnitsTrials_%s_EXP%s_cond%s',runData(iiP).patientName,runData(iiP).EXP,whatToRun.timeStampsType) );
    results = AMNA.getStimulusTriggeredFireRateAllAreas(runData(iiP), fileNameResults, whatToRun);
    PrintActiveFigs(figureFolder)

%     whatToRun.timeStampsType = 'Test';
%     fileNameResults = fullfile(data_p_path, 'populationAnalysis\neuralUnits_trialData', sprintf('neuralUnitsTrials_%s_EXP%s_cond%s',runData(iiP).patientName,runData(iiP).EXP,whatToRun.timeStampsType) );
%     results = AMNA.getStimulusTriggeredFireRateAllAreas(runData(iiP), fileNameResults, whatToRun);
%     PrintActiveFigs(figureFolder)
end

% population - look at activity increment pre-behav point
outputFileFolder = fullfile(data_p_path,'populationAnalysis\microLFP_spectralAnalysis');
outputFigureFolder = fullfile(data_p_path,'populationAnalysis\microLFP_spectralAnalysis\figures');
whatToRun.timeStampsType = 'D';
whatToRun.singlePlots = true;
AMNA.neuralLFP_spectralAnalysis(runData, outputFileFolder, outputFigureFolder, whatToRun);
whatToRun.timeStampsType = 'G';
whatToRun.singlePlots = true;
AMNA.neuralLFP_spectralAnalysis(runData, outputFileFolder, outputFigureFolder, whatToRun);
whatToRun.timeStampsType = 'X';
whatToRun.singlePlots = true;
AMNA.neuralLFP_spectralAnalysis(runData, outputFileFolder, outputFigureFolder, whatToRun);

% cross correlation
% cross-correlation - single-units/multi unit in all brain areas vs
% single-units/multi unit in PROBE area
outputFileFolder = fullfile(data_p_path,'populationAnalysis\neuralUnits_crossCor');
outputFigureFolder = fullfile(data_p_path,'populationAnalysis\neuralUnits_crossCor\figures');
NU_CrossCorOutputFile = fullfile(outputFileFolder,'NUcrossCorrAllPts.mat');
whatTorun.singlePlots = true;
neuralCrossCorrResults = AMNA.neuralCrossCorr(runData, NU_CrossCorOutputFile, outputFigureFolder, whatTorun);

