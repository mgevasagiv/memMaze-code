function   readEdfChannelLocal(edfFilename,edfHdr,ch_id,outputFolder)

data = ft_read_data(edfFilename, 'header',edfHdr,'chanindx',ch_id);
CSC_Sampling_Rate_Hz = 1/edfHdr.Fs;
filename = fullfile(outputFolder, [edfHdr.label{ch_id},'.mat']);
save(filename, 'data','CSC_Sampling_Rate_Hz','edfHdr','-v7.3');

end
