% Pt IR106, MAZE Paradigm
% **************************************************************************** %
%% 1 - Extract NLX MICRO channels, merge files to one joint file named CSC1-80, 
% create averaged version and spike-sort this version 
% Consts contains basic paraemters
dropboxLink
cd(fullfile(dropbox_link,'Code\Nir Lab'))
MatlabPathMaya();

addpath('E:\Dropbox\Code\Useful_code')
addpath('E:\Dropbox\Code\Nir Lab\Work\ONLINE_EEG_TRACKING_analysis')

pt = 'ir106';
exp = 1;
header = getmemMazeExperimentHeader(pt,exp);
preProcessDepthUCD(pt,exp)

poolobj = gcp('nocreate'); % If no pool, do not create new one.
if isempty(poolobj)
    parpool('local',4)
else
    poolsize = poolobj.NumWorkers;
end

dataFolder = 'C:\Users\mgeva\Box\iEEG_MAZE_project\DATA\PATIENT DATA\IR106\Datafiles\Day_2';
edfFilename = fullfile(dataFolder, '2024053015_0002_0001.edf');
edfFilename_part2 = fullfile(dataFolder, '2024053015_0003_0001.edf');

[edfHdr] = ft_read_header(edfFilename);
outputFolder = header.processed_MACRO;

% Save Channel names to Montage-spreadsheet
cnt = 1; chList = {};
for ii = 1:length(edfHdr.label)
    s1 = strfind(edfHdr.label{ii},' ');
    s2 = strfind(edfHdr.label{ii},'-Ref');
    if isempty(s1)||isempty(s2); continue; end    
    chList{cnt} = edfHdr.label{ii}(s1+1:s2-1);
    cnt = cnt + 1;
end
T = array2table(chList','variableNames',{'channels'});
writetable(T, 'E:\MAZE\ptData\macroMontage.csv');

for ch_id = 1:length(edfHdr.label)
    if ~isempty(dir(fullfile(outputFolder,sprintf('%s.mat',edfHdr.label{ch_id})))); continue; end
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

for iif = 2:length(list)
    filename = list(iif).name;
    electName = filename(1:end-4);
    
    % FOR IR-106
    if strcmp(pt,'ir106')
        s1 = strfind(electName,' ');
        s2 = strfind(electName,'-Ref');
        if isempty(s1)||isempty(s2); continue; end
        electName = electName(s1+1:s2-1);
    end

    % parse the elec name
    eName = electName(1:end-1);
    index = str2num(electName(end));
    if strcmp(eName(end),'1')
        eName = electName(1:end-2);
        index = str2num(electName(end-1:end));
    end
        
    chInd = []; cnt = 1;
    while(cnt <= length(MacroMontage))
        a = strcmp(MacroMontage(cnt).Area,eName);
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

% ir106 - data was recorded in two seperate EDF files - need to merge CSCs
% after sub-sampling to 1k Hz.
outputFolder1 = 'E:\MAZE\data_p\pir106\EXP1\Denoised_Downsampled_InMicroVolt\MACRO_PART1';
outputFolder2 = 'E:\MAZE\data_p\pir106\EXP1\Denoised_Downsampled_InMicroVolt\MACRO_PART2';
outputFolder = 'E:\MAZE\data_p\pir106\EXP1\Denoised_Downsampled_InMicroVolt\MACRO';
list1 = dir(fullfile(outputFolder1,'CSC*.mat'));
list2 = dir(fullfile(outputFolder2,'CSC*.mat'));

for ii  = 1:length(list)
    mm = matfile(fullfile(list1(ii).folder,list1(ii).name));
    if ~strcmp(list1(ii).name,list2(ii).name)
        error('files dont match')
    end
    copyfile(fullfile(list1(ii).folder,list1(ii).name),fullfile(outputFolder,list1(ii).name))
    mm_joint = matfile(fullfile(outputFolder,list1(ii).name),'writable',true);
    data1 = mm.data;
    mm1 = matfile(fullfile(list2(ii).folder,list2(ii).name));
    data2 = mm1.data;
    data = [mm.data',mm1.data'];
    mm_joint.data = data;
    mm_joint.merged = true;
end
        