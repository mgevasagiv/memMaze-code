% Pt DA030, MAZE Paradigm
% **************************************************************************** %
%% 1 - Extract NLX MICRO channels, merge files to one joint file named CSC1-80, 
% create averaged version and spike-sort this version 
% Consts contains basic paraemters
dropboxLink
cd(fullfile(dropbox_link,'Code\Nir Lab'))
% MatlabPathMaya();
addpath('E:\Dropbox\Code\Useful_code')

pt = 'da030';
exp = 1;
header = getmemMazeExperimentHeader(pt,exp);
preProcessDepthUCD(pt,exp)

poolobj = gcp('nocreate'); % If no pool, do not create new one.
if isempty(poolobj)
    parpool('local',4)
else
    poolsize = poolobj.NumWorkers;
end

dataFolder = 'C:\Users\mgeva\Box\iEEG_MAZE_project\DATA\PATIENT DATA\DA30\Maze';

outputFolder = header.processed_MACRO;
Ch_name{1} = 'LD';
Ch_name{2} = 'LG';
Ch_name{3} = 'LA';
Ch_name{4} = 'LE';
Ch_name{5} = 'LF';
Ch_name{6} = 'LG';
Ch_name{7} = 'LH';
Ch_name{8} = 'LJ';
Ch_name{9} = 'LK';
Ch_name{10} = 'LL';
Ch_name{11} = 'LN';
Ch_name{12} = 'LM';
Ch_name{13} = 'LP';
Ch_name{14} = 'cu_LG';
Ch_name{15} = 'cu_LD';
Ch_name{16} = 'cu_LG';

parfor ii = 1:length(Ch_name)
    outputFolder = header.processed_MACRO;
    fileList = dir(fullfile(dataFolder,[Ch_name{ii},'*.ncs']));
    PROCESS = [1 0 0 0]; % [Extract files, Merge files, Spike sort]
    try
        getRawCSCSingleFile__NLX(header, dataFolder, outputFolder, fileList, Ch_name{ii}, PROCESS)    
    catch
        continue
    end
end

% Line noise removal is part of CSC extraction for NLX files
% fileList = dir([outputFolder,'\','*.mat']);
% for ii = 1:length(fileList)
%         mm = matfile(fullfile(outputFolder,fileList(ii).name),'Writable',true);
%         fline = 60; % Data recorded in the US
%         samplingRateHz = edfHdr.Fs;
%         denoised = remove_line_noise(mm.data,fline, samplingRateHz);
%         mm.data = denoised;        
%         mm.CSC_Sampling_Rate_Hz = samplingRateHz;
% end


% Rename channels based on macroMontage
mm = matfile(header.macroMontagePath);
MacroMontage = mm.MacroMontage;

list = dir(fullfile(outputFolder,'*.mat'));

for iif = 1:length(list)

    filename = list(iif).name;
    
    electName = filename(1:4);
    if strcmp(electName(1:2),'cu')
        continue
    else
        areaName = electName(1:2);
    end

    if strfind('CSC',areaName)
        continue
    end

    if strcmp(electName(end),'_')
        index = str2num(electName(3));
    else
        index = str2num(electName(3:4));
    end
    
    chInd = []; cnt = 1;
    while(1)
        a = strcmp(MacroMontage(cnt).Name,areaName);
        if a
            chInd = cnt + index - 1;
            if isempty(chInd)
                error('no index')
            end
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
