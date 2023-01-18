classdef AnalyzeMazeNeuralActivity < handle
    
    properties
        samplingRate = 1000;
        
        minSpikeRateToIncludeUnit = 0.1; %Hz ++
        windowSpikeRateAroundStim = 500; %ms ++
        windowSpikeRateForComparison = 200; %ms - for comparing between stim and control
        controlDistForStim = 1000; %ms ++
        firingRateWinSize = 10; %ms ++
        shortTimeRangeAfterStim = 3; %seconds ++
        
    end
    
    methods
        
        function results = getStimulusTriggeredFireRateAllAreas(obj, runData, fileNameResults, whatToRun)
            %
            % The method returns the stimulation triggered spike rate vs control for all the areas of the patient where
            % single unit data exists.
            % 1\21- Updated - The average per channel is an average over the multi-units (and not single units).
            % Multi units are included in the analysis only if the firing rate in them is above minSpikeRateToIncludeUnit.
            % The control per spike is a point controlDistForStim before the stimulations (these are coupled controls).
            % The window around the stimulation (of the returned rate function) is set by windowSpikeRateAroundStim
            % The comparison between stimulation and control (=calculation of p-value) is performed on a window of
            % windowSpikeRateForComparison post stimulation.
            % The rate function is created by using movsum on the events
            % scatter, using a window of avgFigBinning ms.
            %
            % The input runData is a struct in the length of number of patients (for which the analysis is required).
            % In addition it receives the input parameter fileNameResults which includes the file name into which the
            % results will be saved (optional).
            % Each element (=patient) in runData should include the fields:
            % patientName
            % ExpDataFileName – name (including path) of the EXP_DATA for the patient.
            % spikeData – name (including path) of the file which includes the spike data for the patient (i.e. the file
            % which usually has the name <patientName>_spike_timestamps_post_processing.mat).
            %
            % The output struct results includes all the results of the analysis, the plot methods plotStimEffectCortical,
            % plotStimEffectMTL can receive the output as a second results variable and plot the spike rate in the area of
            % the relevant analyzed channel.
            % The output struct is a struct with the length of the number of patients (=the length of runData), where each
            % element includes:
            % patientName
            % areasList – list of the areas for which there is single unit data (according to the spikeData file).
            % nStims – number of stimulations
            % stimTriggeredRates – a cell array in the length of the number of areas, in each element is the average
            % stimulation-triggered spike rate for that area
            % contTriggeredRates - a cell array in the length of the number of areas, in each element is the average
            % control-triggered spike rate for that area
            % channels in that area
            % ratesPerChan - a cell array in the length of the number of areas, where each element is an array of the
            % average spike rates per each channel in that area (over the entire session)
            % nChansInAvg – an array in the length of the number of areas, where each element includes the number of
            % channels included in the average spike rate of that area (i.e. the number of multi units)
            % p – an array in the length of the number of areas, where each element is the p-value comparing between the
            % stimulation triggered and control triggered spike rate for that area. The comparison is between the integral (sum) of the curve after and is calculated using paired t-test.
            
            
            if nargin < 3
                fileNameResults = '';
            end
            
            nPatients = length(runData);
            
            windowSpikeRateForComparison = obj.windowSpikeRateForComparison*obj.samplingRate/1000;
            % indsForSign = obj.windowSpikeRateAroundStim+1:obj.windowSpikeRateAroundStim+windowSpikeRateForComparison;
            indsForSign = obj.windowSpikeRateAroundStim-windowSpikeRateForComparison:obj.windowSpikeRateAroundStim;

            %go over the patients requested
            for iPatient = 1:nPatients
                disp(['Patient ',runData(iPatient).patientName]);
                results(iPatient).patientName = runData(iPatient).patientName;
                
                %load exp data for stimulation timings
                expData = load(runData(iPatient).ExpDataFileName);
                expData = expData.EXP_DATA;
                
                
                spikesRecordingSys = getspikesRecordingSys(str2num(runData(iPatient).patientName(2:end)));
                if strcmp(whatToRun.timeStampsType,'X')
                    stimTimes_forSpikeData = [expData.timestamps_us.X{2}',expData.timestamps_us.X{3}']/1e3;
                    stimTimes_forSpikeData_CONT = [expData.timestamps_us.N{1}',expData.timestamps_us.N{2}',expData.timestamps_us.N{3}']/1e3;
                elseif strcmp(whatToRun.timeStampsType,'D')
                    stimTimes_forSpikeData = [expData.timestamps_us.D{2}',expData.timestamps_us.D{3}']/1e3;
                    stimTimes_forSpikeData_CONT = [expData.timestamps_us.N{1}',expData.timestamps_us.N{2}',expData.timestamps_us.N{3}']/1e3;
                elseif strcmp(whatToRun.timeStampsType,'G')
                    stimTimes_forSpikeData = [expData.timestamps_us.G{2}',expData.timestamps_us.G{3}']/1e3;
                    stimTimes_forSpikeData_CONT = [expData.timestamps_us.N{1}',expData.timestamps_us.N{2}',expData.timestamps_us.N{3}']/1e3;
                end
                stimTimes_forSpikeData_CONT = stimTimes_forSpikeData_CONT(1:end-5);
                
                %check only stim times during NREM
                results(iPatient).stimTimesTest = stimTimes_forSpikeData;
                
                nStims = length(stimTimes_forSpikeData);
                if ~isempty(stimTimes_forSpikeData_CONT)
                    nCont = length(stimTimes_forSpikeData_CONT);
                else
                    nCont = nStims;
                end
                dataDuration = max([stimTimes_forSpikeData,stimTimes_forSpikeData_CONT]) + obj.windowSpikeRateAroundStim*obj.samplingRate*10;
                
                
                %load spike data
                try
                    spikeData = load(runData(iPatient).spikeData);
                catch
                    disp([runData(iPatient).spikeData ' doesn''t exist, continuing']);
                    continue;
                end
                
                if isempty(spikeData.micro_channels_spike_summary.unit_list_xls)
                    continue;
                end
                
                %create a list of all the areas for which there is multi
                %unit data
                allAreas = {spikeData.micro_channels_spike_summary.unit_list_xls.Area};
                areasList = unique(allAreas);
                results(iPatient).areasList = areasList;
                results(iPatient).nStims = nStims;
                nAreas = length(areasList);
                
                channelsList = cell(1,nAreas);
                stimTriggeredRates = cell(1,nAreas);
                contTriggeredRates = cell(1,nAreas);
                ratesPerChan = cell(1,nAreas);
                nChansInAvg = zeros(1,nAreas);
                ps = cell(1,nAreas);
                
                %go over all the areas for the patient
                for iArea = 1:nAreas
                    %find all relevant units
                    unitInds = find(strcmp(allAreas,areasList{iArea}));
                    %find all the channels for that area
                    currChannels = [spikeData.micro_channels_spike_summary.unit_list_xls(unitInds).Channel];
                    channelsList{iArea} = unique(currChannels);
                    %go over channels, merge spikes from that channel,
                    %check if passes rate thresh, and if so find its stim
                    %triggered rate and control triggered rate
                    stimTriggeredRates{iArea} = [];
                    contTriggeredRates{iArea} = [];
                    ratesPerChan{iArea} = zeros(1,length(channelsList{iArea}));
                    stimTriggeredRates{iArea}{1} = zeros(nStims,obj.windowSpikeRateAroundStim*2+1);
                    contTriggeredRates{iArea}{1} = zeros(nCont,obj.windowSpikeRateAroundStim*2+1);
                    
                    nChansInAvg(iArea) = 0;
                    nSUperArea(iArea) = 0;
                    for iChannel = 1:length(channelsList{iArea})
                        currUinds = unitInds(currChannels==channelsList{iArea}(iChannel));
                        nSUperArea(iArea) = nSUperArea(iArea) + length(currUinds);
                        
                        % Get statistics per single unit
                        for iiU = 1:length(currUinds)
                            
                            spikeTimes = sort(cat(1,spikeData.micro_channels_spike_summary.spike_timestamps{currUinds(iiU)}));
                            spikeTimes = spikeTimes(spikeTimes <= dataDuration);
                            [currRateAroundStim, currRateAroundControl] = obj.getStimTriggeredFireRate(spikeTimes, stimTimes_forSpikeData, stimTimes_forSpikeData_CONT);
                            %                             plot(sum(currRateAroundStim),'r'); hold on; plot(sum(currRateAroundControl),'k');
                            %                             legend(sprintf('around %s points',whatToRun.timeStampsType),'control')
                            %                             title(sprintf('Area%d Unit%d',channelsList{iArea}(iChannel),currUinds(iiU)))
                            ratesPerSU{iArea}{iChannel}(iiU) = length(spikeTimes)/(dataDuration/obj.samplingRate); %Hz
                        end
                        
                        %merge spike data to create multi unit
                        spikeTimes = sort(cat(1,spikeData.micro_channels_spike_summary.spike_timestamps{currUinds}));
                        spikeTimes = spikeTimes(spikeTimes <= dataDuration);
                        ratesPerChan{iArea}(iChannel) = length(spikeTimes)/(dataDuration/obj.samplingRate); %Hz
                        
                        %                         check whether the multi unit passes the rate
                        %                         threshold
                        if ratesPerChan{iArea}(iChannel) >= obj.minSpikeRateToIncludeUnit
                            [currRateAroundStim, currRateAroundControl] = obj.getStimTriggeredFireRate(spikeTimes, stimTimes_forSpikeData,stimTimes_forSpikeData_CONT);
                            nChansInAvg(iArea) = nChansInAvg(iArea)+1;
                            stimTriggeredRates{iArea} = stimTriggeredRates{iArea}+currRateAroundStim;
                            contTriggeredRates{iArea} = contTriggeredRates{iArea}+currRateAroundControl;
                        end
                        
                        
                        % Checking response per channel
                        if ratesPerChan{iArea}(iChannel) >= obj.minSpikeRateToIncludeUnit
                            [currRateAroundStim, currRateAroundControl] = obj.getStimTriggeredFireRate(spikeTimes, stimTimes_forSpikeData, stimTimes_forSpikeData_CONT);
                            stimTriggeredRates{iArea}{iChannel}  = currRateAroundStim;
                            contTriggeredRates{iArea}{iChannel} = currRateAroundControl;
                            
                            %(control vs stim)
                            aucStim = sum(stimTriggeredRates{iArea}{iChannel}(:,indsForSign),1)/size(currRateAroundStim,1);
                            aucCont = sum(contTriggeredRates{iArea}{iChannel}(:,indsForSign),1)/size(currRateAroundControl,1);
                            [~,ps{iArea}(iChannel)] = ttest(aucStim,aucCont,'tail','right');

%                             nP = 1000;
%                             disp(sprintf('permutation for area %d',iArea))
%                             [clusters, p_values, t_sums, permutation_distribution ] = permutest( currRateAroundStim', currRateAroundControl', 0, ...
%                                 0.05, nP, true );
%                             clusters_cell{iArea}{iChannel} = clusters;
%                             p_values_cell{iArea}{iChannel} = p_values;
%                             t_sums_cell{iArea}{iChannel} = t_sums;
%                             permutation_distribution_cell{iArea}{iChannel} = permutation_distribution;
                             clear str
%                             if find(p_values < 0.05)
%                                 if sum(currRateAroundStim(clusters{find(p_values < 0.05)})) ~=0 
%                                     str{2} = sprintf('permutation sig. for cluster %d',find(p_values < 0.05));
%                                 end
%                             end
                             str{1} = sprintf('Ch%d MU, ttest = p = %2.2e',channelsList{iArea}(iChannel),ps{iArea}(iChannel));
                            
                            % figure; hold all; for ii = 1:length(clusters); plot(clusters{ii},ii*ones(1,length(clusters{ii})),'r');end
                            figureName = sprintf('%s_A%d_ch%d_stimuli_%s',runData(iPatient).patientName,iArea,iChannel,whatToRun.timeStampsType);
                            figure('name',figureName); shadedErrorBar(1:1001,mean(currRateAroundStim),std(currRateAroundStim)/sqrt(nStims),'lineprops','-r');
                            hold on; shadedErrorBar(1:1001,mean(currRateAroundControl),std(currRateAroundControl)/sqrt(nCont),'lineprops','-k');
                            title(str)
                            plot([500 500],get(gca,'ylim'),'k')
                            axis tight
                            axis([250 750, -inf inf])
                            legend(sprintf('FR around %s points',whatToRun.timeStampsType),'FR control','stimuli')
                            
                        else
                            ps{iArea}(iChannel) = nan;
                        end
                        
                        
                    end
                    
                    % Checking responses per area
                    stimTriggeredRatesAllCh{iArea}= zeros(1,obj.windowSpikeRateAroundStim*2+1);
                    contTriggeredRatesAllCh{iArea}= zeros(1,obj.windowSpikeRateAroundStim*2+1);
                    for iChannel = 1:nChansInAvg
                        stimTriggeredRatesAllCh{iArea} = stimTriggeredRatesAllCh{iArea} + mean(stimTriggeredRates{iArea}{iChannel});
                        contTriggeredRatesAllCh{iArea} = contTriggeredRatesAllCh{iArea} + mean(contTriggeredRates{iArea}{iChannel});
                    end
                    stimTriggeredRatesAllCh{iArea} = stimTriggeredRatesAllCh{iArea}/nChansInAvg;
                    contTriggeredRatesAllCh{iArea} = contTriggeredRatesAllCh{iArea}/nChansInAvg;
                    
                    if nChansInAvg(iArea) > 0
                        %calc pval of area under the curve of the spike rate
                        %(control vs stim)
                        aucStim = stimTriggeredRatesAllCh{iArea}(:,indsForSign);
                        aucCont = contTriggeredRatesAllCh{iArea}(:,indsForSign);
                        [~,psA(iArea)] = ttest(aucStim,aucCont,'tail','right');

                        str{1} = sprintf('aggregated channels, ttest %2.2e',psA(iArea));
                                                
                        figureName = sprintf('allChannelsMU %s A%d stimuliType %s',runData(iPatient).patientName,iArea,whatToRun.timeStampsType);
                        figure('name',figureName)
                        plot(stimTriggeredRatesAllCh{iArea},'r')
                        hold all
                        plot(contTriggeredRatesAllCh{iArea},'k')
                        %shadedErrorBar(1:1001,mean(stimTriggeredRatesAllCh{iArea}),std(stimTriggeredRatesAllCh{iArea})/sqrt(nStims),'lineprops','-r');
                        %hold on; shadedErrorBar(1:1001,mean(contTriggeredRatesAllCh{iArea}),std(contTriggeredRatesAllCh{iArea})/sqrt(nCont),'lineprops','-k');
                        title(str)
                        
                    else
                        ps(iArea) = nan;
                    end
                end
                
                
                
                
                
                results(iPatient).stimTriggeredRates = stimTriggeredRates;
                results(iPatient).contTriggeredRates = contTriggeredRates;
                results(iPatient).channelsList = channelsList;
                results(iPatient).ratesPerChan = ratesPerChan;
                % results(iPatient).ratesPerSU = ratesPerSU;
                results(iPatient).nSUperArea = nSUperArea;
                results(iPatient).nChansInAvg = nChansInAvg;
                results(iPatient).p = ps;
                %                 results(iPatient).clusters_cell = clusters_cell;
                %                 results(iPatient).p_values_cell = p_values_cell;
                %                 results(iPatient).permutation_distribution_cell = permutation_distribution_cell;
                %                 results(iPatient).t_sums_cell = t_sums_cell;
                
                PLOT = 0;
                if PLOT
                    if ~isempty(results(iPatient).areasList)
                        iArea = 1;
                        figureName = sprintf('%s A%d stimuli %s',runData(iPatient).patientName,iArea,iChannel,whatToRun.timeStampsType);
                        figure('name',figureName)
                        for iiC = 1:length(results(iPatient).stimTriggeredRates{iArea})
                            resStim = results(iPatient).stimTriggeredRates{iArea}{iiC};
                            resStimControl = results(iPatient).contTriggeredRates{iArea}{iiC};
                            if ~isempty(resStim)
                                subplot(2,4,iiC);
                                shadedErrorBar([-obj.windowSpikeRateAroundStim:obj.windowSpikeRateAroundStim]/obj.samplingRate, mean(resStim), std(resStim)/sqrt(nStims),'lineprops','-m');
                                hold all;
                                shadedErrorBar([-obj.windowSpikeRateAroundStim:obj.windowSpikeRateAroundStim]/obj.samplingRate, mean(resStimControl), std(resStimControl)/sqrt(nStims),'lineprops','-k');
                                % legend({'Stim','Control'});
                                xlabel('Time (sec)');
                                ylabel('Spike rate');
                                chansInAvg = 1; % results(iPatient).channelsList{iArea}(results(iPatient).ratesPerChan{iArea}>=obj.minSpikeRateToIncludeUnit);
                                Parea = results(iPatient).p{iArea};
                                title({['Area: ',results(iPatient).areasList{iArea}],['Chans in avg: ',num2str(chansInAvg)],['pval (stim>cont) - ',num2str(Parea(iiC))]});
                            end
                        end
                    end
                end
                
                
                
                
            end
            
            if ~isempty(fileNameResults)
                save(fileNameResults,'results');
            end
            
        end
        
        
        function [rateAroundStim, rateAroundControl] = getStimTriggeredFireRate(obj, spikeTimes, stimTimes, controlTimes)
            % Receives as input spike times (an array) and stimulation times (an array) and returns:
            % rateAroundStim – a matrix where each row is the rate function around a stimulation.
            % rateAroundControl – a matrix where each row is the rate function around a control.
            %
            %
            dataDuration = max(stimTimes) + obj.shortTimeRangeAfterStim*obj.samplingRate;
            spikeTimes = spikeTimes(spikeTimes<=dataDuration);
            
            nStim = length(stimTimes);
            
            %convert spike times to rate function using movsum
            spikes = zeros(1, floor(dataDuration));
            spikes(round(spikeTimes)) = 1;
            spikeRateSession = movsum(spikes,obj.firingRateWinSize,'Endpoints','fill')/(obj.firingRateWinSize/1000);
            
            
            %go over stimulations and save the rate function around the
            %stimulation and around control
            
            %             % Option 1 - constant delay
            %             for iStim =1:nStim
            %
            %                 currStimTime = stimTimes(iStim);
            %                 currControlTime = currStimTime-obj.controlDistForStim;
            %
            %                 rateAroundStim(iStim,:) = spikeRateSession(currStimTime-obj.windowSpikeRateAroundStim:currStimTime+obj.windowSpikeRateAroundStim);
            %                 rateAroundControl(iStim,:) = spikeRateSession(currControlTime-obj.windowSpikeRateAroundStim:currControlTime+obj.windowSpikeRateAroundStim);
            %
            %             end
            
            if isempty(controlTimes)
                rateAroundStim = zeros(nStim,obj.windowSpikeRateAroundStim*2+1);
                rateAroundControl = zeros(nStim,obj.windowSpikeRateAroundStim*2+1);
                
                % Option 2 - random delay from stim time
                randShiftMax = 500;
                for iStim =1:nStim
                    
                    currStimTime = stimTimes(iStim);
                    randShift = rand(1);
                    currControlTime = currStimTime - (randShift*randShiftMax + obj.controlDistForStim);
                    
                    rateAroundStim(iStim,:) = spikeRateSession(floor(currStimTime-obj.windowSpikeRateAroundStim):floor(currStimTime+obj.windowSpikeRateAroundStim));
                    rateAroundControl(iStim,:) = spikeRateSession(floor(currControlTime-obj.windowSpikeRateAroundStim):floor(currControlTime+obj.windowSpikeRateAroundStim));
                    
                end
            else
                
                rateAroundStim = zeros(nStim,obj.windowSpikeRateAroundStim*2+1);
                rateAroundControl = zeros(length(controlTimes),obj.windowSpikeRateAroundStim*2+1);
                
                
                for iStim =1:nStim
                    currStimTime = stimTimes(iStim);
                    rateAroundStim(iStim,:) = spikeRateSession(floor(currStimTime-obj.windowSpikeRateAroundStim):floor(currStimTime+obj.windowSpikeRateAroundStim));
                end
                for iStim =1:length(controlTimes)-2
                    currStimTime = controlTimes(iStim);
                    rateAroundControl(iStim,:) = spikeRateSession(floor(currStimTime-obj.windowSpikeRateAroundStim):floor(currStimTime+obj.windowSpikeRateAroundStim));
                end
            end
            
            
        end
        
        
    end
    
end

