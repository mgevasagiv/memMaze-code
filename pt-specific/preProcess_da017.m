% Pt DA017, MAZE Paradigm
% **************************************************************************** %
%% 1 - Extract NLX MICRO channels, merge files to one joint file named CSC1-80, 
% create averaged version and spike-sort this version 
% Consts contains basic paraemters
dropboxLink
cd(fullfile(dropbox_link,'Code\Nir Lab'))
MatlabPathMaya();

pt = 'pDA017';
exp = 1;
header = getmemMazeExperimentHeader(pt,exp);

poolobj = gcp('nocreate'); % If no pool, do not create new one.
if isempty(poolobj)
    parpool('local',4)
else
    poolsize = poolobj.NumWorkers;
end

dataFolder = 'E:\MY_BACKUP\UCD_raw_datasets\PATIENT DATA\DA17\Maze_Maya\';
outputFolder = header.processed_MACRO;
Ch_name{1} = 'SMG';
Ch_name{2} = 'PTRI';
Ch_name{3} = 'LES';
Ch_name{4} = 'IPOS';
Ch_name{5} = 'IPRE';
Ch_name{6} = 'PINS';
Ch_name{7} = 'cu_ptri';
Ch_name{8} = 'MPOS';
Ch_name{9} = 'MPRE';
Ch_name{10} = 'AINS';
Ch_name{11} = 'ACIN';
Ch_name{12} = 'cu_acin';

parfor ii = 1:length(Ch_name)
    outputFolder = header.processed_MACRO;
    fileList = dir(fullfile(dataFolder,[Ch_name{ii},'*.ncs']));
    PROCESS = [1 0 0 0]; % [Extract files, Merge files, Spike sort]
    getRawCSCSingleFile__NLX(header, dataFolder, outputFolder, fileList, Ch_name{ii}, PROCESS)    
end

%%
header.excel_sheet = 'E:\MAZE\ptData\da17.xlsx';
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

%%
header.spikes_sampling_Rate_Hz = 32e3;
% Spike sorting was ran with SR updated to *?kHz* in wave_clus
spikeSort_maze(pt,exp)
