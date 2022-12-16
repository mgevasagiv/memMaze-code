% Pt DA018, MAZE Paradigm
% **************************************************************************** %
%% 1 - Extract NLX MICRO channels, merge files to one joint file named CSC1-80, 
% create averaged version and spike-sort this version 
% Consts contains basic paraemters
dropboxLink
cd(fullfile(dropbox_link,'Code\Nir Lab'))
MatlabPathMaya();

pt = 'pda018';
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


%% Read (manually set) montage from XLS - based on blue lab notebook --> mat file
macroMontageSheet = header.macro_montage_sheet;
p_in.numeric_fields = {'Channel','Electrode_number'}; % numeric fields in excel
cell_rows = 1:EXP_DATA.ELECTRODES*8+8;
MacroMontage= [];

if strfind(header.id,'p') % manual montage for UCLA patients
    MacroMontage = read_excel_sheet(header.excel_sheet,macroMontageSheet,cell_rows,p_in.numeric_fields,p_in);
    
    % get rid of empty entries
    ii = length(MacroMontage);
    while(1)
        if isempty(MacroMontage(ii).Area)
            MacroMontage(ii) = [];
            disp(sprintf('removed entry %d',ii))
        else
            break
        end
        ii = ii - 1;
    end
    if isempty(dir(header.root_data))
        mkdir(header.root_data)
    end
    save(header.macroMontagePath,'MacroMontage')
end

%%