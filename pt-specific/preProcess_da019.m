% Pt DA017, MAZE Paradigm
% **************************************************************************** %
%% 1 - Extract NLX MICRO channels, merge files to one joint file named CSC1-80, 
% create averaged version and spike-sort this version 
% Consts contains basic paraemters
dropboxLink
cd(fullfile(dropbox_link,'Code\Nir Lab'))
MatlabPathMaya();

pt = 'da019';
exp = 1;
header = getmemMazeExperimentHeader(pt,exp);
if isempty(header.processedDataPath)
    mkdir(header.processedDataPath)
end
preProcessDepthUCD(pt,exp)

poolobj = gcp('nocreate'); % If no pool, do not create new one.
if isempty(poolobj)
    parpool('local',4)
else
    poolsize = poolobj.NumWorkers;
end

datasetFilename = fullfile(header.processedDataPath,sprintf('p%s_EXP%d_dataset.mat',pt,exp));
load(datasetFilename)


dataFolder = 'C:\Users\mgeva\Box\MAZE_project\DATA\PATIENT DATA\DA19\MayaMaze_1';
outputFolder = header.processed_MACRO;
Ch_name{1} = 'cu_lam';
Ch_name{2} = 'cu_ldpi';
Ch_name{3} = 'LAM';
Ch_name{4} = 'LAMFG';
Ch_name{5} = 'LDPCI';
Ch_name{6} = 'LDPI';
Ch_name{7} = 'LDPRE';
Ch_name{8} = 'LMMFG';
Ch_name{9} = 'LOFC';
Ch_name{10} = 'LPA';
Ch_name{11} = 'LPL';
Ch_name{12} = 'LPMFG';
Ch_name{13} = 'Photodiode';

outputFolder = header.processed_MACRO;
parfor ii = 1:length(Ch_name)
    fileList = dir(fullfile(dataFolder,[Ch_name{ii},'*.ncs']));
    PROCESS = [1 0 0 0]; % [Extract files, Merge files, Spike sort]
    getRawCSCSingleFile__NLX(header, dataFolder, outputFolder, fileList, Ch_name{ii}, PROCESS)    
end

% Rename channels based on macroMontage
mm = matfile(header.macroMontagePath);
MacroMontage = mm.MacroMontage;

list = dir(fullfile(outputFolder,'*.mat'));

for iif = 11:length(list)
    filename = list(iif).name;
    electName = filename(1:end-9);
    
    % parse the elec name
    areaName = electName(1:end-1);
    index = str2num(electName(end));
    if strcmp(areaName(end),'1')
        areaName = electName(1:end-2);
        index = str2num(electName(end-1:end));
    end
        
    chInd = []; cnt = 1;
    while(cnt < length(MacroMontage))
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


%% 
%%
header.excel_sheet = 'E:\MAZE\ptData\da19.xlsx';
header.micro_montage_sheet = 'microMontage';
EXP_DATA.MICRO_ELECTRODES = 2;

microMontageSheet = header.micro_montage_sheet;
p_in.numeric_fields = {'Channel'}; % numeric fields in excel
cell_rows = 1:EXP_DATA.MICRO_ELECTRODES*8+8;
MicroMontage= [];

if strfind(header.id,'p') % manual montage for UCLA patients
    MicroMontage = read_excel_sheet(header.excel_sheet,microMontageSheet,cell_rows,p_in.numeric_fields,p_in);
    
    % get rid of empty entries
    ii = length(MicroMontage);
    while(1)
        if isempty(MicroMontage(ii).Area)
            MicroMontage(ii) = [];
            disp(sprintf('removed entry %d',ii))
        else
            break
        end
        ii = ii - 1;
    end
    if isempty(dir(header.root_data))
        mkdir(header.root_data)
    end
    save(header.microMontagePath,'MicroMontage')
end

datasetFilename = fullfile(header.processedDataPath,sprintf('%s_EXP%d_dataset.mat',header.id,header.experimentNum));
save(datasetFilename,'header')

%% Feb 9th, 2023
header.spikes_sampling_Rate_Hz = 32e3;
% Spike sorting was ran with SR updated to *?kHz* in wave_clus
spikeSort_maze(pt,exp)


