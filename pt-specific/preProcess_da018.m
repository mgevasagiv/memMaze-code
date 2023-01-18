% Pt DA018, MAZE Paradigm
% **************************************************************************** %
%% 1 - Extract NLX MICRO channels, merge files to one joint file named CSC1-80, 
% create averaged version and spike-sort this version 
% Consts contains basic paraemters
dropboxLink
cd(fullfile(dropbox_link,'Code\Nir Lab'))
MatlabPathMaya();

pt = 'da018';
exp = 1;
header = getmemMazeExperimentHeader(pt,exp);

poolobj = gcp('nocreate'); % If no pool, do not create new one.
if isempty(poolobj)
    parpool('local',4)
else
    poolsize = poolobj.NumWorkers;
end

dataFolder = 'E:\MY_BACKUP\UCD_raw_datasets\PATIENT DATA\DA18\Maze_Maya\';



outputFolder = header.processed_MACRO;
Ch_name{1} = 'RASL';
Ch_name{2} = 'RAIL';
Ch_name{3} = 'RPIL';
Ch_name{4} = 'RPSL';
Ch_name{5} = 'RMC';
Ch_name{6} = 'RAC';
Ch_name{7} = 'LAC';
Ch_name{8} = 'LAMC';
Ch_name{9} = 'LPMC';
Ch_name{10} = 'PHOTODIODE';
Ch_name{11} = 'AUDIO';


edfFilename = 'E:\MY_BACKUP\UCD_raw_datasets\PATIENT DATA\DA18\EPHYS_Natus\Maya_Maze\DA018.edf';
[edfHdr] = ft_read_header(edfFilename);
outputFolder = header.processed_MACRO;

parfor ch_id = 1:length(edfHdr.label)
    if ~contains(edfHdr.label{ch_id},'eeg','IgnoreCase',true); continue; end
    readEdfChannelLocal(edfFilename,edfHdr,ch_id,outputFolder)  
end

fileList = dir([outputFolder,'\','*.mat']);
for ii = 1:length(fileList)
        mm = matfile(fullfile(outputFolder,fileList(ii).name),'Writable',true);
        fline = 60; % Data recorded in the US
        samplingRateHz = 1e3;
        denoised = remove_line_noise(mm.data,fline, samplingRateHz);
        mm.data = denoised;        
        mm.CSC_Sampling_Rate_Hz = samplingRateHz;
end


% Rename channels based on macroMontage

mm = matfile(header.macroMontagePath);
MacroMontage = mm.MacroMontage;

list = dir(fullfile(outputFolder,'*.mat'));

for iif = 1:length(list)
    filename = list(iif).name;
    electName = filename(5:end-8);
    
    % parse the elec name
    areaName = electName(1:end-1);
    index = str2num(electName(end));
    if strcmp(areaName(end),'1')
        areaName = electName(1:end-2);
        index = str2num(electName(end-1:end));
    end
        
    chInd = []; cnt = 1;
    while(1)
        a = strcmp(MacroMontage(cnt).Area,areaName);
        if a
            chInd = cnt + index - 1;
            fprintf('rename %s to %s\n',fullfile(outputFolder,filename),fullfile(outputFolder,sprintf('CSC%d.mat',chInd)))
            disp('')
            movefile(fullfile(outputFolder,filename),fullfile(outputFolder,sprintf('CSC%d.mat',chInd)));
            cnt = 0; clear chInd
            break
        end
        cnt = cnt + 1;
    end
end