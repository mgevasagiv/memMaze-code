% Pt DA023, MAZE Paradigm
% **************************************************************************** %
%% 1 - Extract NLX MICRO channels, merge files to one joint file named CSC1-80, 
% create averaged version and spike-sort this version 
% Consts contains basic paraemters
dropboxLink
cd(fullfile(dropbox_link,'Code\Nir Lab'))
MatlabPathMaya();

pt = 'IR103';
exp = 1;
header = getmemMazeExperimentHeader(pt,exp);
preProcessDepthUCD(pt,exp)

poolobj = gcp('nocreate'); % If no pool, do not create new one.
if isempty(poolobj)
    parpool('local',4)
else
    poolsize = poolobj.NumWorkers;
end

dataFolder = 'C:\Users\mgeva\Box\MAZE_project\DATA\PATIENT DATA\IR103\2023060110';
besaFilename = fullfile(dataFolder,'2023060110_0001.besa');

dataFolder = 'C:\Users\mgeva\Box\MAZE_project\DATA\PATIENT DATA\IR103\2023060111';
besaFilename = fullfile(dataFolder,'2023060111_0001.besa');

outputFolder = header.processed_MACRO;

hdr = ft_read_header(besaFilename);
cfg = [];
cfg.dataset = besaFilename;
data_struct = ft_preprocessing(cfg);
TRIAL = cell2mat(data_struct.trial);

for ch_id = 1:length(hdr.label)
    filename = hdr.label{ch_id};
    data = TRIAL(ch_id,:);
    CSC_Sampling_Rate_Hz = hdr.Fs;
    save(fullfile(outputFolder,filename),'data','hdr','CSC_Sampling_Rate_Hz','-v7.3')
end

fileList = dir([outputFolder,'\','*.mat']);
for ii = 1:length(fileList)
        mm = matfile(fullfile(outputFolder,fileList(ii).name),'Writable',true);
        fline = 60; % Data recorded in the US
        samplingRateHz = hdr.Fs;
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
