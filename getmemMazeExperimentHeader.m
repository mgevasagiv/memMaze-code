% Based on Amit M. Code - adapted by MayaG - July 2016
function header = getmemMazeExperimentHeader(patientId, experimentNum)

header.id = patientId;

[~, hostname]= system('hostname');
if strfind(hostname, 'cns-cdkngw2')
    root_data_general = 'E:\MAZE\Data_p\MACRO_MONTAGE';
    root_datap_general = 'E:\MAZE\Data_p';
elseif strfind(hostname, 'Laptop')
    root_data_general = 'C:\Users\MAZE\Data\';
    root_datap_general =  'C:\Users\MAZE\Data_p\';
else
    disp('ERROR: Computer name not recognized')
    return
end

%if no experiment num than this is the default is the first experiment
if nargin==2
    experimentNum = 1;
end

header.experimentNum = experimentNum;
root_data = fullfile(root_data_general, sprintf('%s',patientId),sprintf('EXP%d',experimentNum));
root_processed = fullfile(root_datap_general, sprintf('%s',patientId),sprintf('EXP%d',experimentNum));
header.root_data = root_data;
header.isAutoMontage = false;
header.microMontagePath = fullfile(root_data,'microMontage.mat');
header.macroMontagePath = fullfile(root_data,'MacroMontage.mat');
header.processedDataPath = fullfile(root_processed);
header.figuresDataPath =  fullfile(root_processed,'figures');
header.spikesDataPath = fullfile(root_processed,'spikeSorting');
header.processed_MACRO = fullfile(root_processed,'Denoised_Downsampled_InMicroVolt','MACRO');
header.processed_MICRO = fullfile(root_processed,'Denoised_Downsampled_InMicroVolt','MICRO');
if isempty(header.figuresDataPath); mkdir(header.figuresDataPath); end
if isempty(header.processed_MACRO); mkdir(header.processed_MACRO); end
if isempty(header.processed_MICRO); mkdir(header.processed_MICRO); end
if isempty(header.spikesDataPath); mkdir(header.spikesDataPath); end

%%---------------------------------------------------------------------------------------------%%

if (strcmp(lower(patientId),'da017')) % Feb 25th 2016 at Ichilov, pre-defined MACRO Stimulation during Nap - Firas and Yuval
    if isAwake
    else
        if (experimentNum==1)
            header.StimExp = 1; % default
            header.StimType = 'MACRO';
            
            header.nsFilePath{Index.MACRO} = '';
            header.nsFormat{Index.MACRO} = 'ns3';
            header.channelsVec{Index.MACRO} = 1:99;
            header.denoisedChannelsPath{Index.MACRO} = '';
            header.dataIndexInNsFile{Index.MACRO} = 2; % index of biggest data chunk
            
            header.nsFilePath{Index.MICRO} = '';
            header.nsFormat{Index.MICRO} = 'ns5';
            header.channelsVec{Index.MICRO} = 1:88;
            header.denoisedChannelsPath{Index.MICRO} = '';
            header.dataIndexInNsFile{Index.MICRO} = 2; % index of biggest data chunk
            
            header.nevFilePath = '';
            
            header.SpikeSortingXLS = ''; % all identified units
            
        end
    end
 


else 
    
    header.isValid = false;
    
end % patients

end
