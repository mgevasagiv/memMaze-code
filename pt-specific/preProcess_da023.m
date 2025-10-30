% Pt DA023, MAZE Paradigm
% **************************************************************************** %
%% 1 - Extract NLX MICRO channels, merge files to one joint file named CSC1-80, 
% create averaged version and spike-sort this version 
% Consts contains basic paraemters
dropboxLink
cd(fullfile(dropbox_link,'Code\Nir Lab'))
MatlabPathMaya();

pt = 'da023';
exp = 1;
header = getmemMazeExperimentHeader(pt,exp);
preProcessDepthUCD(pt,exp)

poolobj = gcp('nocreate'); % If no pool, do not create new one.
if isempty(poolobj)
    parpool('local',4)
else
    poolsize = poolobj.NumWorkers;
end

dataFolder = 'E:\MY_BACKUP\UCD_raw_datasets\PATIENT DATA\';

edfFilename = 'E:\MY_BACKUP\UCD_raw_datasets\PATIENT DATA\DA23\DA23_Maze.EDF';
[edfHdr] = ft_read_header(edfFilename);
outputFolder = header.processed_MACRO;

parfor ch_id = 1:length(edfHdr.label)
    % if ~contains(edfHdr.label{ch_id},'eeg','IgnoreCase',true); continue; end
    readEdfChannelLocal(edfFilename,edfHdr,ch_id,outputFolder)  
end

fileList = dir([outputFolder,'\','*.mat']);
for ii = 1:length(fileList)
        mm = matfile(fullfile(outputFolder,fileList(ii).name),'Writable',true);
        fline = 60; % Data recorded in the US
        samplingRateHz = edfHdr.Fs;
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
    electName = filename(1:end-4);
    
    if(strfind(filename,'-'))
        % parse the elec name
        doubleDigit = false;
        areaName = electName(1:end-2);
        if strcmp(areaName(end),'-')
            areaName = electName(1:end-3);
            doubleDigit = true;
        end
        if ~doubleDigit
            index = str2num(electName(end));
        else
            index = str2num(electName(end-1:end));
        end
    else
        doubleDigit = false;
        areaName = electName(1:end-1);
        if strcmp(areaName(end),'1')
            areaName = electName(1:end-2);
            doubleDigit = true;
        end
        if ~doubleDigit
            index = str2num(electName(end));
        else
            index = str2num(electName(end-1:end));
        end
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

list = dir(fullfile(outputFolder,'CSC*.mat'));
LFP_SamplingRate_Hz = 1000;
for ii = 1:length(list)
    mm = matfile(fullfile(outputFolder, list(ii).name),'Writable',true);
    if mm.CSC_Sampling_Rate_Hz == LFP_SamplingRate_Hz; continue; end
    data = accurateResampling(mm.data, mm.CSC_Sampling_Rate_Hz, LFP_SamplingRate_Hz);
    mm.data = data;
    mm.CSC_Sampling_Rate_Hz = LFP_SamplingRate_Hz;
    clear mm
end
