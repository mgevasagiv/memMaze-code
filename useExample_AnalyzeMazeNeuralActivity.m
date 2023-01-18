patients = {'pda017'};
expNames = {'EXP1'};

data_p_path = 'E:\MAZE\Data_p\';
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
    
    results_folder_top = 'E:\Data_p\ClosedLoopDataset\';

end


AMNA = AnalyzeMazeNeuralActivity;
fileNameResults = [];
whatToRun.timeStampsType = 'D';
results = AMNA.getStimulusTriggeredFireRateAllAreas(runData, fileNameResults, whatToRun);
PrintActiveFigs('E:\MAZE\data_p\pda017\EXP1\Figures')

whatToRun.timeStampsType = 'X';
results = AMNA.getStimulusTriggeredFireRateAllAreas(runData, fileNameResults, whatToRun);
PrintActiveFigs('E:\MAZE\data_p\pda017\EXP1\Figures')

whatToRun.timeStampsType = 'G';
results = AMNA.getStimulusTriggeredFireRateAllAreas(runData, fileNameResults, whatToRun);
PrintActiveFigs('E:\MAZE\data_p\pda017\EXP1\Figures')