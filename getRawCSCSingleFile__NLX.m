% Maya GS 
function getRawCSCSingleFile__NLX(header, dataFolder, outputFolder, files, Ch_name, PROCESS, target_folder_DOWNSAMPLE, details)

if PROCESS(1)
    
    for ii = 1:length(files)
        
        filename_CSC_EEG_in = fullfile(dataFolder,files(ii).name);
        filename = fullfile(outputFolder,sprintf('%s.mat',files(ii).name(1:end-4)));
        
        a = dir(filename);
        if (~isempty(a))
            disp(sprintf('file %s already extracted',files(ii).name))
            continue
        end
        disp(sprintf('extracting %s',filename_CSC_EEG_in))
        
        FieldSelection = [1 1 1 1 1] ; % Will read all the variables from the file
        ExtractionMode = 2 ; % read only relevant session according to given timestamps
        ExtractionModeArray = [1 600000];
        
        a = dir(filename_CSC_EEG_in);
        if( a.bytes == 16384 )
            disp(sprintf('file %s is empty (%d bytes)',files(ii).name,files(ii).bytes))
            continue
        end
        [Timestamps, ChanNum, SampleFrequency, NumValidSamples, Samples, NlxHeader] = ...
            Nlx2MatCSC( filename_CSC_EEG_in, FieldSelection, 1, ExtractionMode, ExtractionModeArray ) ;
        if ~NumValidSamples
            disp(sprintf('file %s is empty (%d bytes)',files(ii).name,files(ii).bytes))
            continue
        end
        if length(unique(ChanNum)) > 1 || length(unique(SampleFrequency)) > 1
            warning('Something wrong with the NLX acquisition')
        end
        % check how many samples were dropped
        MISS = sum(Samples(:) == 0);
        timeLost = MISS/SampleFrequency(1);
        percLost = 100*MISS/sum(NumValidSamples);
        disp(sprintf('%2.2f sec dropped in rec, %2.3f%%%',timeLost,percLost))
        CSC_SamplePeriod_microsec = round( 10^6/SampleFrequency(1) );
        % CSC_SamplePeriod_microsec =  10^6/SampleFrequency(1);
        CSC_Sampling_Rate_Hz = SampleFrequency(1) ;
        Samples_reshaped = zeros(1, prod(size(Samples)) );
        Timestamps_filledIn = zeros(1, prod(size(Samples)) );
        for ii_DataBlock = 1:size(Samples,2), % Loop over the 512-point blocks of data
            idx_data = (1:size(Samples,1)) + (ii_DataBlock-1)*size(Samples,1); % Indexes where to put the data
            dataChunk = NaN(1,size(Samples,1));
            if NumValidSamples(ii_DataBlock) > 0
                dataChunk(1:NumValidSamples(ii_DataBlock)) = Samples(1:NumValidSamples(ii_DataBlock),ii_DataBlock)';
            end
            Samples_reshaped( idx_data ) = dataChunk;
            Timestamps_filledIn( idx_data ) = Timestamps(ii_DataBlock) + (1:size(Samples,1))*CSC_SamplePeriod_microsec;
        end
        % This is tricky - when NumValidSamples < 512 the rest of the data
        % in the matching column is junk (a copy of previous column), here
        % I remove these excess data-points -
        Samples = Samples_reshaped;
        Timestamps = Timestamps_filledIn; % Timestamps in microsec
        
        Samples(isnan(Samples_reshaped)) = 0;
        Timestamps(isnan(Samples_reshaped)) = 0;
        
        if ( length(Samples)-sum(isnan(Samples_reshaped)) ~= sum(NumValidSamples) )
            error('extraction failed')
        end
        
        % Aug 2017, Oct 2017 patients had a dropped-data error, so I'm
        % going to remove all zero-valued samples BFR applying any
        % additional adjustments to the signal
        rmvIdx = (Samples == 0);
        Samples(rmvIdx) = [];
        Timestamps(rmvIdx) = [];
        raw_data_interp = interp1(Timestamps, Samples, Timestamps_filledIn);
        raw_data_interp(isnan(raw_data_interp)) = 0; % assuming that only the edges can be NaN at this stage
        disp(sprintf('%2.2f sec (%2.2f %%) interpolated',(length(Timestamps_filledIn) - length(Timestamps))/CSC_Sampling_Rate_Hz,...
                                        100*(length(Timestamps_filledIn) - length(Timestamps))/length(Timestamps_filledIn)))
        Samples = raw_data_interp - nanmean(raw_data_interp); % The data, with Mean Removed
        Timestamps = Timestamps_filledIn;
        
        SampleFrequency = SampleFrequency(1);
        clear  Timestamps_filledIn  Samples_reshaped  ; % Clear some large and unnecessary variables, to avoid "OUT OF MEMORY" problems
        % Convert the data ("Samples") to Microvolts and filter out (band-stop) the 60-Hz / 120 Hz noises:
        Samples_microvolts = Samples * str2num( NlxHeader{17}(14:end) ) * 10^6 ; % Convert to Microvolts
        % Do not decimate uchannels data
        % decimate_factor = 1/(CSC_SamplePeriod_microsec/round( 10^6/Consts.DOWNSAMPLED_FREQUENCY )) ;
        % Samples_microvolts_decimated = decimate( Samples_microvolts, decimate_factor ); % Decimate
        
        fline = 60; % Data recorded in the US
        denoised = remove_line_noise(Samples_microvolts,fline, SampleFrequency(1));
        data = denoised;
        
        if (length(Samples) ~= length(Timestamps))
            error('Extraction failed')
        end
        save(filename, 'data','Timestamps','CSC_Sampling_Rate_Hz','NlxHeader','-v7.3')
        
    end
    
end


end