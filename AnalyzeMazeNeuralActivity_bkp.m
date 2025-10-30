classdef AnalyzeMazeNeuralActivity < handle
    
    properties
        samplingRate = 1000;
        
        minSpikeRateToIncludeUnit = 0.1; %Hz ++
        windowSpikeRateAroundStim = 500; %ms ++
        windowSpikeRateForComparison = 200; %ms - for comparing between stim and control
        controlDistForStim = 1000; %ms ++
        firingRateWinSize = 10; %ms ++
        shortTimeRangeAfterStim = 3; %seconds ++
        
        windowLfpForComparison = 200; % ms
        % auto-correlation vector size
        corrBuffer = 250; %msec
        
        crossCorfsz = 8;
    end
    
    methods
        
        function MACROLFP_collectPerMazeNeuralMarkers(obj, runData, outputFileFolder, outputFigureFolder, whatToRun)
            global data_p_path
            global BEHAV_TRIALS_FOLDER
            
            if nargin < 3
                fileNameResults = '';
            end
            
            nPatients = length(runData);
            
            
            %go over the patients requested
            ich_cnt = 0;
            for iPatient = 1:nPatients
                
                clear AREAS
                
                disp(['Patient ',runData(iPatient).patientName]);
                results.patientName = runData(iPatient).patientName;
                results.EXP = runData(iPatient).EXP;
                
                
                filename =  fullfile(BEHAV_TRIALS_FOLDER, [results.patientName(2:end),'_',results.EXP,'_data_structure.mat']);
                behavfile = dir(filename);
                try
                    mm = matfile(fullfile(BEHAV_TRIALS_FOLDER, behavfile(1).name));
                catch
                    disp('behav file missing')
                    return;
                end
                
                try
                    eval(['data_structure = mm.da',results.patientName(4:end),'_data_structure;'])
                    test_data_structure = mm.test_data_structure;
                catch
                    eval(['data_structure = mm.ir',results.patientName(4:end),'_data_structure;'])
                    test_data_structure = mm.test_data_structure;
                end
                
                % behavior markers
                filename =  fullfile(BEHAV_TRIALS_FOLDER, [results.patientName(2:end),'_',results.EXP,'behav_data_structure.mat']);
                behavfile = dir(filename);
                mBehav = matfile(fullfile(BEHAV_TRIALS_FOLDER, behavfile(1).name));
                behav_data_structure = mBehav.behav_data_structure;
                mazeCompletionTable = mBehav.mazeCompletionTable;
                
                list = 0:23;
                clear sorted_rep id_rep
                bType = whatToRun.timeStampsType;
                mazeIds =  eval(['behav_data_structure.',bType,'_rep1_maze_id']);
                missingMazes = find(~ismember(list,unique(mazeIds,'stable')));

%                 if strcmp(results.patientName(2:end),  'da018') & strcmp(results.EXP,  '2')
%                     list = 0:23;
%                     missingMazes = find(~ismember(list,(behav_data_structure.G_rep1_maze_id)));
%                     % missingMazes = [missingMazes 11];
%                 end
                mazeCompletionTable.solve_ms_rep1(missingMazes) = NaN;
              
                [id_rep{1}, sorted_rep{1}] = sort(unique(mazeIds,'stable'));
                [id_rep{2}, sorted_rep{2}] = sort(unique(mazeIds,'stable'));
                [id_rep{3}, sorted_rep{3}] = sort(unique(mazeIds,'stable'));
                if ~isempty(missingMazes)
                    for ii = 1:3
                        rmv = ismember(list(missingMazes),id_rep{ii});
                        id_rep{ii}(missingMazes(rmv)) = [];
                        sorted_rep{ii}(missingMazes(rmv)) = [];
                    end
                end
                
                
                obj.samplingRate = data_structure.sampling_rate;
                Nelectrodes = size(data_structure.D_rep1,3);
                aucContrastAcrossRep1Rep2 = nan(1, Nelectrodes);
                aucContrastAcrossRep1Rep3 = nan(1, Nelectrodes);
                try
                    disp(runData(iPatient).macroMontageFileName)
                    load(runData(iPatient).macroMontageFileName)
                catch
                    error('macro montage file not found')
                end
                
                for ii = 1:Nelectrodes
                    contact_id = data_structure.channels_included(ii);
                    AREAS{ii} = MacroMontage(contact_id).Area;
                end
                
                if strcmp(whatToRun.timeStampsType,'X')
                    trials1 = data_structure.X_rep1;
                    cont_trials1 = data_structure.N_rep1;
                    trials2 = data_structure.X_rep2;
                    cont_trials2 = data_structure.N_rep2;
                    trials3 = data_structure.X_rep3;
                    cont_trials3 = data_structure.N_rep3;
                    
                    idx2 = data_structure.X_rep2_good_maze_trials;
                    idx3 = data_structure.X_rep3_good_maze_trials;
                    trialsMem = [data_structure.X_rep2(idx2,:,:); data_structure.X_rep3(idx3,:,:)];
                    idx2 = data_structure.N_rep2_good_maze_trials;
                    idx3 = data_structure.N_rep3_good_maze_trials;
                    cont_trialsMem =  [data_structure.N_rep2(idx2,:,:); data_structure.N_rep3(idx3,:,:)];
                    
                elseif strcmp(whatToRun.timeStampsType,'G')
                    trials1 = data_structure.G_rep1;
                    cont_trials1 = data_structure.N_rep1;
                    trials2 = data_structure.G_rep2;
                    cont_trials2 = data_structure.N_rep2;
                    trials3 = data_structure.G_rep3;
                    cont_trials3 = data_structure.N_rep3;
                    
                    idx2 = data_structure.G_rep2_good_maze_trials;
                    idx3 = data_structure.G_rep3_good_maze_trials;
                    trialsMem = [data_structure.G_rep2(idx2,:,:); data_structure.G_rep3(idx3,:,:)];
                    idx2 = data_structure.N_rep2_good_maze_trials;
                    idx3 = data_structure.N_rep3_good_maze_trials;
                    cont_trialsMem =  [data_structure.N_rep2(idx2,:,:); data_structure.N_rep3(idx3,:,:)];
                    
                elseif strcmp(whatToRun.timeStampsType,'D')
                    trials1 = data_structure.D_rep1;
                    cont_trials1 = data_structure.N_rep1;
                    trials2 = data_structure.D_rep2;
                    cont_trials2 = data_structure.N_rep2;
                    trials3 = data_structure.D_rep3;
                    cont_trials3 = data_structure.N_rep3;
                    
                    idx2 = data_structure.D_rep2_good_maze_trials;
                    idx3 = data_structure.D_rep3_good_maze_trials;
                    trialsMem = [data_structure.D_rep2(idx2,:,:); data_structure.D_rep3(idx3,:,:)];
                    idx2 = data_structure.N_rep2_good_maze_trials;
                    idx3 = data_structure.N_rep3_good_maze_trials;
                    cont_trialsMem =  [data_structure.N_rep2(idx2,:,:); data_structure.N_rep3(idx3,:,:)];
                end
                
                
                windowLfpForComparisonBfrAftr_ms = floor(obj.windowLfpForComparison*obj.samplingRate/1000);
                midPoint = floor(size(trials1,2)/2);
                indsForSign = [midPoint-windowLfpForComparisonBfrAftr_ms:midPoint+windowLfpForComparisonBfrAftr_ms];
                
                for iiC = 1:Nelectrodes
                    disp(sprintf('%d/%d electrodes',iiC,Nelectrodes))
                    
                    
                    % aucStim = abs(sum(trials1(:,indsForSign,iiC)'))/length(indsForSign);
                    % aucCont = abs(sum(cont_trials1(:,indsForSign,iiC)'))/length(indsForSign);
                    
                    pkMax = max(trials1(:,indsForSign,iiC)');
                    pkMin = min(trials1(:,indsForSign,iiC)');
                    aucStim = pkMax - pkMin;
                    pkMax = max(cont_trials1(:,indsForSign,iiC)');
                    pkMin = min(cont_trials1(:,indsForSign,iiC)');
                    aucCont = pkMax - pkMin;
                    
                    try
                        [psA,~] = ranksum(aucStim,aucCont);
                    catch
                        pSA = 1;
                        warning('rank sum was not ran')
                    end
                    %aucStim = mean(abs(sum(trials1(:,indsForSign,iiC)')))/length(indsForSign);
                    %aucCont = mean(abs(sum(cont_trials1(:,indsForSign,iiC)')))/length(indsForSign);
                    aucStim = mean(aucStim);
                    aucCont = mean(aucCont);
                    
                    psA1(iiC) = psA;
                    aucStimRep1(iiC) = aucStim;
                    aucContRep1(iiC) = aucCont;
                    aucContrast1(iiC) = (aucStim-aucCont)./(aucStim+aucCont);
                    aucRatio1(iiC,:) = 100*(aucStim-aucCont)./max(aucStim,aucCont);
                    
                    % REP2
                    aucStim = abs(sum(trials2(:,indsForSign,iiC)'))/length(indsForSign);
                    aucCont = abs(sum(cont_trials2(:,indsForSign,iiC)'))/length(indsForSign);
                    try
                        [psA,~] = ranksum(aucStim,aucCont);
                    catch
                        pSA = 1;
                        warning('rank sum was not ran')
                    end
                    aucStim = mean(abs(sum(trials2(:,indsForSign,iiC)')))/length(indsForSign);
                    aucCont = mean(abs(sum(cont_trials2(:,indsForSign,iiC)')))/length(indsForSign);
                    psA2(iiC) = psA;
                    aucStimRep2(iiC) = aucStim;
                    aucContRep2(iiC) = aucCont;
                    aucContrast2(iiC) = (aucStim-aucCont)./(aucStim+aucCont);
                    aucRatio2(iiC,:) = 100*(aucStim-aucCont)./max(aucStim,aucCont);
                    
                    % REP3
                    aucStim = abs(sum(trials3(:,indsForSign,iiC)'))/length(indsForSign);
                    aucCont = abs(sum(cont_trials3(:,indsForSign,iiC)'))/length(indsForSign);
                    try
                        [psA,~] = ranksum(aucStim,aucCont);
                    catch
                        pSA = 1;
                        warning('rank sum was not ran')
                    end
                    aucStim = mean(abs(sum(trials3(:,indsForSign,iiC)')))/length(indsForSign);
                    aucCont = mean(abs(sum(cont_trials3(:,indsForSign,iiC)')))/length(indsForSign);
                    psA3(iiC) = psA;
                    aucStimRep3(iiC) = aucStim;
                    aucContRep3(iiC) = aucCont;
                    aucContrast3(iiC) = (aucStim-aucCont)./(aucStim+aucCont);
                    aucRatio3(iiC,:) = 100*(aucStim-aucCont)./max(aucStim,aucCont);
                    
                    % MEM
                    aucStim = (sum(trialsMem(:,indsForSign,iiC)'))/length(indsForSign);
                    aucCont = (sum(cont_trialsMem(:,indsForSign,iiC)'))/length(indsForSign);
                    try
                        [psA,~] = ranksum(aucStim,aucCont);
                    catch
                        pSA = 1;
                        warning('rank sum was not ran')
                    end
                    aucStim = mean(sum(trialsMem(:,indsForSign,iiC)'))/length(indsForSign);
                    aucCont = mean(sum(cont_trialsMem(:,indsForSign,iiC)'))/length(indsForSign);
                    psA4(iiC) = psA;
                    aucStimRep4(iiC) = aucStim;
                    aucContRep4(iiC) = aucCont;
                    aucContrast4(iiC) = (aucStim-aucCont)./(aucStim+aucCont);
                    aucRatio4(iiC,:) = 100*(aucStim-aucCont)./max(aucStim,aucCont);
                    
                    
                    % Looking at how different channels activity change
                    % across trials
                    if strcmp(whatToRun.timeStampsType,'D')
                        aucStim1 = mean(abs(sum(trials1(:,indsForSign,iiC)')))/length(indsForSign);
                        aucStim2 = mean(abs(sum(trials2(:,indsForSign,iiC)')))/length(indsForSign);
                        aucStim3 = mean(abs(sum(trials3(:,indsForSign,iiC)')))/length(indsForSign);
                        aucRatioAcrossRep1Rep2(iiC) = mean(100*(aucStim1(:)-aucStim2(:))./max(aucStim1',aucStim2'));
                        aucRatioAcrossRep1Rep3(iiC) = mean(100*(aucStim1(:)-aucStim3(:))./max(aucStim1',aucStim3'));
                        aucContrastAcrossRep1Rep2(iiC) = mean((aucStim1(:)-aucStim2(:))./(aucStim1(:)+aucStim2(:)));
                        aucContrastAcrossRep1Rep3(iiC) = mean((aucStim1(:)-aucStim3(:))./(aucStim1(:)+aucStim3(:)));
                    end
                    
                    if sum(strcmp(whatToRun.timeStampsType,{'G','X'}))
                        aucStim1 = sum(abs(trials1(sorted_rep{1},indsForSign,iiC)')/length(indsForSign));
                        aucStim2 = sum(abs(trials2(sorted_rep{2},indsForSign,iiC)')/length(indsForSign));
                        aucStim3 = sum(abs(trials3(sorted_rep{3},indsForSign,iiC)')/length(indsForSign));
                        aucRatioAcrossRep1Rep2(iiC) = nanmean(100*(aucStim1(:)-aucStim2(:))./max(aucStim1',aucStim2'));
                        aucRatioAcrossRep1Rep3(iiC) = nanmean(100*(aucStim1(:)-aucStim3(:))./max(aucStim1',aucStim3'));
                        aucContrastAcrossRep1Rep2(iiC) = nanmean((aucStim1(:)-aucStim2(:))./(aucStim1(:)+aucStim2(:)));
                        aucContrastAcrossRep1Rep3(iiC) = nanmean((aucStim1(:)-aucStim3(:))./(aucStim1(:)+aucStim3(:)));
                        
                        % Correlate with behavior
                        solve_ms_rep1 = mazeCompletionTable.solve_ms_rep1;
                        solve_ms_rep1(missingMazes) = [];
                        solve_ms_rep1(isnan(solve_ms_rep1)) = 0;
                        
                        solve_ms_rep2 = mazeCompletionTable.solve_ms_rep2;
                        solve_ms_rep2(missingMazes) = [];
                        solve_ms_rep2(isnan(solve_ms_rep2)) = 0;
                        
                        solve_ms_rep3 = mazeCompletionTable.solve_ms_rep3;
                        solve_ms_rep3(missingMazes) = [];
                        solve_ms_rep3(isnan(solve_ms_rep3)) = 0;
                        
                        RT_diff_rep1_rep2 = [solve_ms_rep2 - solve_ms_rep1]/1e6;
                        RT_diff_rep2_rep3 = [solve_ms_rep3 - solve_ms_rep2]/1e6;
                        aucContrastAcrossRep1Rep2_allMazes(iiC,:) = ((aucStim1(:)-aucStim2(:))./(aucStim1(:)+aucStim2(:)));
                        aucContrastAcrossRep2Rep3_allMazes(iiC,:) = ((aucStim2(:)-aucStim3(:))./(aucStim2(:)+aucStim3(:)));
                        try
                            [R,p] = corrcoef(RT_diff_rep1_rep2,aucContrastAcrossRep1Rep2_allMazes(iiC,:));
                            P_corr(iiC,1) = p(1,2);
                            [R,p] = corrcoef(RT_diff_rep2_rep3,aucContrastAcrossRep2Rep3_allMazes(iiC,:));
                            P_corr(iiC,2) = p(1,2);
                        catch
                            disp('corr failed')
                        end
                        
                    end
                    
                end
                
                if sum(strcmp(whatToRun.timeStampsType,{'G','X'}))
                    ind1 = find(RT_diff_rep1_rep2 > 0);
                    M = aucContrastAcrossRep1Rep2_allMazes(:,ind1);
                    ind2 = find(RT_diff_rep1_rep2 < 0);
                    N = aucContrastAcrossRep1Rep2_allMazes(:,ind2);
                    [p_improve, h] = ranksum(M(:),N(:));
                    disp(sprintf('%s - P - neural-behav between improved and degraded routes %2.2e', whatToRun.timeStampsType, p_improve))
                end
                
                for ii = 1:length(AREAS)
                    area = classifyArea_MAZE(AREAS{ii});
                    if area.isHip
                        hipCh(ii) = 1;
                    else
                        hipCh(ii) = 0;
                    end
                    if area.isMTL
                        MTLCh(ii) = 1;
                    else
                        MTLCh(ii) = 0;
                    end
                    if area.isFrontal
                        frontCh(ii) = 1;
                    else
                        frontCh(ii) = 0;
                    end
                end
                
                ptID = ones(1,length(psA1))*iPatient;
                TablePerPt = table(ptID(:), data_structure.channels_included(:),AREAS', frontCh(:), MTLCh(:), hipCh(:), ...
                    aucContrast1(:), aucStimRep1(:) , psA1(:),...
                    aucContrast2(:), aucStimRep2(:) , psA2(:),...
                    aucContrast3(:), aucStimRep3(:) , psA3(:),...
                    aucContrast4(:), aucStimRep4(:) , psA4(:),...
                    aucRatio1(:), aucRatio2(:),aucRatio3(:),aucRatio4(:),....
                    aucRatioAcrossRep1Rep2(:),aucRatioAcrossRep1Rep3(:),...
                    aucContrastAcrossRep1Rep2(:),aucContrastAcrossRep1Rep3(:),...
                    'VariableNames',{'ptID','chan','AREAS','isFrontal','isMTL','isHip',...
                    'aucCI1','trialAuc1','pSA1',...
                    'aucCI2','trialAuc2','pSA2',...
                    'aucCI3','trialAuc3','pSA3',...
                    'aucCIMem','trialAucMem','pSAMem',...
                    'aucRatio1','aucRatio2','aucRatio3','aucRatioMem',...
                    'aucRatioAcrossRep1Rep2','aucRatioAcrossRep1Rep3',...
                    'aucContrastAcrossRep1Rep2','aucContrastAcrossRep1Rep3'});
                
                %                 TablePerPt_Ch = table(allCh_aucContrast1);
                
                if iPatient == 1
                    TABLE_ALL = TablePerPt;
                else
                    
                    TABLE_ALL = [TABLE_ALL; TablePerPt];
                end
                
                if sum(strcmp(whatToRun.timeStampsType,{'G','X'}))
                    filename = sprintf('%s_%s_behavNeural_%s',runData(iPatient).patientName,runData(iPatient).EXP, whatToRun.timeStampsType);
                    fileNameResults = fullfile(outputFileFolder, filename);
                    channelsIncluded = data_structure.channels_included(:);
                    save(fileNameResults,'p_improve','aucContrastAcrossRep1Rep2_allMazes','aucContrastAcrossRep2Rep3_allMazes','mazeCompletionTable','AREAS','channelsIncluded','P_corr'...
                        ,'frontCh','MTLCh','hipCh','missingMazes','RT_diff_rep2_rep3','RT_diff_rep1_rep2')
                end
                                    
                clear AREAS aucStimRep1 aucStimRep2 aucContrast1 aucContrast2 psA1 psA2 hipCh MTLCh frontCh
                clear aucContrast3 aucStimRep3 psA3
                clear aucRatio1 aucRatio2 aucRatio3 aucRatio4
                clear aucContrast4 aucStimRep4 psA4
                clear aucRatioAcrossRep1Rep2 aucRatioAcrossRep1Rep3 aucContrastAcrossRep1Rep2 aucContrastAcrossRep1Rep3
                clear aucContrastAcrossRep1Rep2_allMazes aucContrastAcrossRep2Rep3_allMazes
            end
            
            
            
            if ~isempty(outputFileFolder)
                filename = sprintf('allPts_stimuli_%s_MACRO_aucTable',whatToRun.timeStampsType);
                fileNameResults = fullfile(outputFileFolder, filename);
                save(fileNameResults,'TABLE_ALL');
            end
            
        end
        
        
        
        function neuralLFP_spectralAnalysis(obj, runData, outputFileFolder, outputFigureFolder, whatToRun)
            
            global data_p_path
            
            if nargin < 3
                fileNameResults = '';
            end
            
            nPatients = length(runData);
            
            
            %go over the patients requested
            for iPatient = 1:nPatients
                disp(['Patient ',runData(iPatient).patientName]);
                results(iPatient).patientName = runData(iPatient).patientName;
                microChannels_spectralAnalysis = runData(iPatient).microChannels_spectralAnalysis;
                
                %load exp data for stimulation timings
                expData = load(runData(iPatient).ExpDataFileName);
                expData = expData.EXP_DATA;
                
                for iiC = 1:length(microChannels_spectralAnalysis);
                    
                    try
                        channel = microChannels_spectralAnalysis(iiC);
                        mm = matfile(fullfile(data_p_path, 'pda017\EXP1\spikeSorting',sprintf('csc%d.mat',channel)));
                        data = mm.data;
                        obj.samplingRate = mm.CSC_Sampling_Rate_Hz;
                        results(iPatient).spectralAnalysis(iiC).channelID = channel;
                    catch
                        disp('lfp data not found')
                    end
                    windowLfpForComparisonBfrAftr_ms = obj.windowLfpForComparison*obj.samplingRate/1000;
                    % indsForSign = obj.windowSpikeRateAroundStim-windowSpikeRateForComparison:obj.windowSpikeRateAroundStim;
                    
                    if strcmp(whatToRun.timeStampsType,'X')
                        stimTimes = [expData.timestamps_us.X{2}',expData.timestamps_us.X{3}']/1e3;
                        stimTimes_CONT = [expData.timestamps_us.N{1}',expData.timestamps_us.N{2}',expData.timestamps_us.N{3}']/1e3;
                        
                    elseif strcmp(whatToRun.timeStampsType,'D')
                        stimTimes = [expData.timestamps_us.D{2}',expData.timestamps_us.D{3}']/1e3;
                        stimTimes_CONT = [expData.timestamps_us.N{1}',expData.timestamps_us.N{2}',expData.timestamps_us.N{3}']/1e3;
                        
                    elseif strcmp(whatToRun.timeStampsType,'G')
                        stimTimes = [expData.timestamps_us.G{2}',expData.timestamps_us.G{3}']/1e3;
                        stimTimes_CONT = [expData.timestamps_us.N{1}',expData.timestamps_us.N{2}',expData.timestamps_us.N{3}']/1e3;
                    end
                    stimTimes_forSpikeData_CONT = stimTimes_CONT(1:end-5);
                    results(iPatient).stimTimesTest = stimTimes;
                    results(iPatient).stimTimes_forSpikeData_CONT = stimTimes_CONT;
                    
                    nStims = length(stimTimes);
                    if ~isempty(stimTimes_CONT)
                        nCont = length(stimTimes_CONT);
                    else
                        nCont = nStims;
                    end
                    dataDuration = max([stimTimes,stimTimes_CONT]) + obj.windowSpikeRateAroundStim*obj.samplingRate*10;
                    
                    clear trials cont_trials
                    for ii_s = 1:length(stimTimes)
                        trials(ii_s,:) = data(stimTimes(ii_s) - windowLfpForComparisonBfrAftr_ms:stimTimes(ii_s) + windowLfpForComparisonBfrAftr_ms);
                    end
                    for ii_s = 1:length(stimTimes_CONT)
                        cont_trials(ii_s,:) = data(stimTimes_CONT(ii_s) - windowLfpForComparisonBfrAftr_ms:stimTimes_CONT(ii_s) + windowLfpForComparisonBfrAftr_ms);
                    end
                    
                    time_halfbandwidth = 2.5;
                    dps_seq = dpss(size(trials,1),time_halfbandwidth);
                    params.Fs = obj.samplingRate;
                    params.trialave  = true;
                    params.fpass = [0 200];
                    p = 0.05;
                    params.err  = [2 p]; % Jackknife error bars
                    % estimate the spectrum of averaged single trials based on micro local field potential
                    % trials expected as times * trials
                    [estSpectrum,freq,confidenceBands] = mtspectrumc( trials', params );
                    [estSpectrum_cont,freq,confidenceBands_cont] = mtspectrumc( cont_trials', params );
                    results(iPatient).spectralAnalysis(iiC).estSpectrum = estSpectrum;
                    results(iPatient).spectralAnalysis(iiC).estSpectrum_cont = estSpectrum_cont;
                    results(iPatient).spectralAnalysis(iiC).freq = freq;
                    results(iPatient).spectralAnalysis(iiC).confidenceBands = confidenceBands;
                    results(iPatient).spectralAnalysis(iiC).spectrumCalcParams = params;
                    
                    % Add moving window to estimate the spectrum of averaged single trials based on micro local field potential
                    % Note units have to be consistent. Thus, if movingwin is in seconds, Fs
                    % has to be in Hz. see chronux.m for more information.
                    clear params
                    params.Fs = obj.samplingRate; % Hz
                    params.fpass = [0 200];
                    params.trialave  = true;
                    movingwin = [0.1 0.01] ; % sec
                    params.tapers=[3 5]; % very smoothed (default was [3 5]
                    [estSpectrum_movingwin,t_movingwin,freq_movingwin] = mtspecgramc( trials', movingwin, params );
                    [estSpectrum_movingwin_cont,t_movingwin,freq_movingwin] = mtspecgramc( cont_trials', movingwin, params );
                    
                    results(iPatient).spectralAnalysis(iiC).estSpectrum_movingwin = estSpectrum_movingwin;
                    results(iPatient).spectralAnalysis(iiC).estSpectrum_movingwin_cont = estSpectrum_movingwin_cont;
                    results(iPatient).spectralAnalysis(iiC).t_movingwin = t_movingwin;
                    results(iPatient).spectralAnalysis(iiC).freq_movingwin = freq_movingwin;
                    results(iPatient).spectralAnalysis(iiC).spectrumCalcParams.params = params;
                    results(iPatient).spectralAnalysis(iiC).spectrumCalcParams.movingwin = movingwin;
                    
                    if whatToRun.singlePlots
                        
                        figureName = sprintf('%s_ch%d_stimuli_%s_spectralAnalysis',runData(iPatient).patientName,channel,whatToRun.timeStampsType);
                        newA4figure(figureName)
                        subplot(3,2,1:2); hold all
                        plot(freq,10*log10(estSpectrum),'b')
                        plot(freq,10*log10(estSpectrum_cont),'k')
                        legend(whatToRun.timeStampsType,'N')
                        
                        subplot(3,2,3)
                        CLIM = [prctile(10*log10(estSpectrum_movingwin(:)),15) prctile(10*log10(estSpectrum_movingwin(:)),85) ];
                        h1 = imagesc(t_movingwin,freq_movingwin, 10*log10(estSpectrum_movingwin),CLIM);
                        axis([0.1 0.5 0 200])
                        set(gca,'xlim',[0.1 0.5], 'XTick',[0.1 0.3 0.5],'XTickLabel',{'-0.2', '0', '0.2'})
                        hold all
                        plot(0.3*ones(1,2),get(gca,'ylim'),'k')
                        colorbar
                        title(sprintf('trials - %s',whatToRun.timeStampsType))
                        axis ij
                        ylabel('Frequency (Hz)')
                        xlabel('Time (sec)')
                        
                        subplot(3,2,4)
                        h2 = imagesc(t_movingwin,freq_movingwin, 10*log10(estSpectrum_movingwin_cont),CLIM);
                        title(sprintf('control trials - %s','N'))
                        set(gca,'xlim',[0.1 0.5], 'XTick',[0.1 0.3 0.5],'XTickLabel',{'-0.2', '0', '0.2'})
                        colorbar
                        hold all
                        plot(0.3*ones(1,2),get(gca,'ylim'),'k')
                        xlabel('Time (sec)')
                        axis ij
                        
                        PrintActiveFigs(outputFigureFolder)
                        
                    end
                end
                
            end % patients
            
            if ~isempty(outputFileFolder)
                filename = sprintf('%s_stimuli_%s_spectralAnalysis',runData(iPatient).patientName,whatToRun.timeStampsType);
                fileNameResults = fullfile(outputFileFolder, filename);
                save(fileNameResults,'results');
            end
            
        end % function
        
        
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
                
                if strcmp(whatToRun.timeStampsType,'X')
                    stimTimes_forSpikeData = [expData.timestamps_us.X{2}',expData.timestamps_us.X{3}']/1e3;
                    stimTimes_forSpikeData_CONT = [expData.timestamps_us.N{1}',expData.timestamps_us.N{2}',expData.timestamps_us.N{3}']/1e3;
                    
                elseif strcmp(whatToRun.timeStampsType,'D')
                    stimTimes_forSpikeData = [expData.timestamps_us.D{2}',expData.timestamps_us.D{3}']/1e3;
                    % stimTimes_forSpikeData = [expData.timestamps_us.D{3}']/1e3;
                    stimTimes_forSpikeData_CONT = [expData.timestamps_us.N{1}',expData.timestamps_us.N{2}',expData.timestamps_us.N{3}']/1e3;
                    
                elseif strcmp(whatToRun.timeStampsType,'G')
                    stimTimes_forSpikeData = [expData.timestamps_us.G{2}',expData.timestamps_us.G{3}']/1e3;
                    stimTimes_forSpikeData_CONT = [expData.timestamps_us.N{1}',expData.timestamps_us.N{2}',expData.timestamps_us.N{3}']/1e3;
                    
                    
                elseif strcmp(whatToRun.timeStampsType,'Gm')
                    stimTimes_forSpikeData = [expData.timestamps_us.G{2}',expData.timestamps_us.G{3}']/1e3;
                    stimTimes_forSpikeData_CONT = [expData.timestamps_us.G{1}']/1e3;
                    
                elseif strcmp(whatToRun.timeStampsType,'Test') %%DBG
                    trials = test_data_structure.trials_correct_choices(:,1:1000,:);
                    trials(isnan(trials)) = 0;
                    
                    idx2 = data_structure.N_rep2_good_maze_trials;
                    idx3 = data_structure.N_rep3_good_maze_trials;
                    cont_trials =  [data_structure.N_rep2(idx2,:,:); data_structure.N_rep3(idx3,:,:)];
                    
                    stimTimes_forSpikeData = [expData.timestamps_us.G{2}',expData.timestamps_us.G{3}']/1e3;
                    stimTimes_forSpikeData_CONT = [expData.timestamps_us.G{1}']/1e3;
                    
                end
                stimTimes_forSpikeData_CONT = stimTimes_forSpikeData_CONT(1:end-5);
                results(iPatient).stimTimesTest = stimTimes_forSpikeData;
                results(iPatient).stimTimes_forSpikeData_CONT = stimTimes_forSpikeData_CONT;
                
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
                    stimTriggeredRates_SU{iArea}{1}{1} = zeros(nCont,obj.windowSpikeRateAroundStim*2+1);
                    
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
                            stimTriggeredRates_SU{iArea}{iChannel}{iiU} =  currRateAroundStim;
                        end
                        
                        %merge spike data to create multi unit
                        spikeTimes = sort(cat(1,spikeData.micro_channels_spike_summary.spike_timestamps{currUinds}));
                        spikeTimes = spikeTimes(spikeTimes <= dataDuration);
                        ratesPerChan{iArea}(iChannel) = length(spikeTimes)/(dataDuration/obj.samplingRate); %Hz
                        
                        %                         check whether the multi unit passes the rate
                        %                         threshold
                        %                         if ratesPerChan{iArea}(iChannel) >= obj.minSpikeRateToIncludeUnit
                        %                             [currRateAroundStim, currRateAroundControl] = obj.getStimTriggeredFireRate(spikeTimes, stimTimes_forSpikeData,stimTimes_forSpikeData_CONT);
                        %                             nChansInAvg(iArea) = nChansInAvg(iArea)+1;
                        %                             stimTriggeredRates{iArea} = stimTriggeredRates{iArea}+currRateAroundStim;
                        %                             contTriggeredRates{iArea} = contTriggeredRates{iArea}+currRateAroundControl;
                        %                         end
                        
                        
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
                            str{1} = sprintf('Ch%d MU [-200,0]ms: ttest = p = %2.2e',channelsList{iArea}(iChannel),ps{iArea}(iChannel));
                            
                            
                            figureName = sprintf('%s_ch%d_stimuli_%s',runData(iPatient).patientName,channelsList{iArea}(iChannel),whatToRun.timeStampsType);
                            newA4figure(figureName)
                            subplot(3,1,1); hold all
                            h = imagesc(currRateAroundStim);
                            plot([500 500],get(gca,'ylim'),'w')
                            colormap parula
                            set(gca,'xlim',[250 750], 'XTick',[250 500 750],'XTickLabel',{'-0.25', '0', '0.25'},...
                                'ylim',[1 size(currRateAroundStim,1)])
                            ylabel('trials')
                            axis ij
                            title_str{1} =  sprintf('%s, %s, ch%d',runData(iPatient).patientName,areasList{iArea},channelsList{iArea}(iChannel));
                            title_str{2} = sprintf('rate around %s points',whatToRun.timeStampsType);
                            title(title_str)
                            
                            subplot(3,1,2); hold all
                            h = imagesc(currRateAroundControl);
                            plot([500 500],get(gca,'ylim'),'w')
                            set(gca,'xlim',[250 750], 'XTick',[250 500 750],'XTickLabel',{'-0.25', '0', '0.25'},...
                                'ylim',[1 size(currRateAroundControl,1)])
                            ylabel('trials')
                            axis ij
                            title('rate around control points')
                            
                            subplot(3,1,3)
                            p1 = shadedErrorBar(1:1001,mean(currRateAroundStim),std(currRateAroundStim)/sqrt(nStims),'lineprops','-r');
                            hold on; p2 = shadedErrorBar(1:1001,mean(currRateAroundControl),std(currRateAroundControl)/sqrt(nCont),'lineprops','-k');
                            
                            title(str)
                            plot([500 500],get(gca,'ylim'),'k')
                            axis tight
                            YLIM = get(gca,'ylim');
                            % plot(indsForSign([1,end]),ones(1,2)*YLIM(2)-0.05*diff(YLIM),'c','linewidth',3)
                            starLoc = max(mean(currRateAroundStim))*1.2;
                            
                            idx2test = 250:750;
                            for iiT = 1:length(idx2test)
                                a = currRateAroundStim(:,idx2test(iiT));
                                b = currRateAroundControl(:,idx2test(iiT));
                                
                                [T df] = ttest2_cell( { a b } ,'inhomogenous');
                                signif(iiT) = 2*tcdf(-abs(T), df(1));
                                if signif(iiT) < 0.05
                                    plot(idx2test(iiT),starLoc,'b*');
                                end
                            end
                            
                            axis([250 750, -inf inf])
                            set(gca,'XTick',[250 500 750],'XTickLabel',{'-0.25', '0', '0.25'})
                            legend([p1.mainLine,p2.mainLine], {sprintf('FR around %s points',whatToRun.timeStampsType),'FR control'})
                            
                        else
                            ps{iArea}(iChannel) = nan;
                        end
                        
                        
                    end
                    
                    % Checking responses per area
                    stimTriggeredRatesAllCh{iArea}= zeros(1,obj.windowSpikeRateAroundStim*2+1);
                    contTriggeredRatesAllCh{iArea}= zeros(1,obj.windowSpikeRateAroundStim*2+1);
                    for iChannel = 1:nChansInAvg(iArea)
                        stimTriggeredRatesAllCh{iArea} = stimTriggeredRatesAllCh{iArea} + mean(stimTriggeredRates{iArea}{iChannel});
                        contTriggeredRatesAllCh{iArea} = contTriggeredRatesAllCh{iArea} + mean(contTriggeredRates{iArea}{iChannel});
                    end
                    stimTriggeredRatesAllCh{iArea} = stimTriggeredRatesAllCh{iArea}/nChansInAvg(iArea);
                    contTriggeredRatesAllCh{iArea} = contTriggeredRatesAllCh{iArea}/nChansInAvg(iArea);
                    
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
                        psA(iArea) = nan;
                    end
                end
                
                
                
                
                results(iPatient).stimTriggeredRates_SU = stimTriggeredRates_SU;
                results(iPatient).stimTriggeredRates = stimTriggeredRates;
                results(iPatient).contTriggeredRates = contTriggeredRates;
                results(iPatient).channelsList = channelsList;
                results(iPatient).ratesPerChan = ratesPerChan;
                results(iPatient).ratesPerSU = ratesPerSU;
                results(iPatient).nSUperArea = nSUperArea;
                results(iPatient).nChansInAvg = nChansInAvg;
                results(iPatient).p = ps;
                
                
                PLOT = 1;
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
                    if  stimTimes(iStim) > length(spikeRateSession)
                        break
                    end
                    currStimTime = stimTimes(iStim);
                    rateAroundStim(iStim,:) = spikeRateSession(floor(currStimTime-obj.windowSpikeRateAroundStim):floor(currStimTime+obj.windowSpikeRateAroundStim));
                end
                for iStim =1:length(controlTimes)-2
                    if  controlTimes(iStim) > length(spikeRateSession)
                        break
                    end
                    currStimTime = controlTimes(iStim);
                    rateAroundControl(iStim,:) = spikeRateSession(floor(currStimTime-obj.windowSpikeRateAroundStim):floor(currStimTime+obj.windowSpikeRateAroundStim));
                end
            end
            
            
        end
        
        
        
        function results = neuralCrossCorr(obj, runData, fileNameResults, folderToSaveFigures, whatTorun)
            
            if nargin < 3
                fileNameResults = '';
            end
            
            global data_p_path
            nPatients = length(runData);
            
            for iPatient = 1:nPatients
                disp(['Patient ',runData(iPatient).patientName]);
                
                results(iPatient).patientName = runData(iPatient).patientName;
                results(iPatient).EXP = runData(iPatient).EXP;
                pt = runData(iPatient).patientName;
                
                %load exp data for stimulation timings
                expData = load(runData(iPatient).ExpDataFileName);
                expData = expData.EXP_DATA;
                
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
                nAreas = length(areasList);
                
                channelsList = cell(1,nAreas);
                
                
                %go over all the areas for the patient
                for iArea = 1:nAreas
                    %find all relevant units
                    unitInds = find(strcmp(allAreas,areasList{iArea}));
                    %find all the channels for that area
                    currChannels = [spikeData.micro_channels_spike_summary.unit_list_xls(unitInds).Channel];
                    channelsList{iArea} = unique(currChannels);
                    
                    nChansInAvg(iArea) = 0;
                    conditions = {'D','X','G','N'};
                    
                    for ii_c = 1:length(conditions)
                        
                        if ii_c ~= 4
                            mm = matfile(fullfile('E:\MAZE\data_p\populationAnalysis\neuralUnits_trialData\',...
                                sprintf('neuralUnitsTrials_%s_EXP%s_cond%s.mat',runData(iPatient).patientName,runData(iPatient).EXP,conditions{ii_c})));
                            trialDef = mm.results;
                            stimTimes = trialDef.stimTimesTest;
                        else
                            stimTimes = trialDef.stimTimes_forSpikeData_CONT;
                        end
                        nstims = length(stimTimes);
                        
                        
                        for iChannel = 1:length(channelsList{iArea})
                            
                            currUinds = unitInds(currChannels==channelsList{iArea}(iChannel));
                            
                            % For each single unit - look at changes in phase
                            % locking
                            for ii_su = 1:length(currUinds)
                                
                                su_ind = currUinds(ii_su);
                                spikeTimes = spikeData.micro_channels_spike_summary.spike_timestamps{su_ind}(:)';
                                
                                spikeInfo.spikeTimes_ms =  spikeData.micro_channels_spike_summary.spike_timestamps{ii_su};
                                spikeInfo.spike_shapes_mean =  spikeData.micro_channels_spike_summary.spike_shapes_mean{ii_su};
                                spikeInfo.spike_shapes_std =  spikeData.micro_channels_spike_summary.spike_shapes_std{ii_su};
                                spikeInfo.su_ind = su_ind;
                                results(iPatient).spikeCELL{iArea, iChannel,ii_su}.spikeInfo = spikeInfo;
                                
                                events_ms_e = [];
                                for ii_e = 1:nstims
                                    t0 = stimTimes(ii_e);
                                    events_ms_e = [events_ms_e, spikeTimes(spikeTimes > t0-obj.corrBuffer & spikeTimes < t0+obj.corrBuffer)];
                                end
                                results(iPatient).spikeCELL{iArea, iChannel,ii_su}.spikeTimes_ms{ii_c} = events_ms_e; % spikes in entire epoch
                            end % unit
                            
                        end % channel
                        
                    end % condition
                    
                    
                end % area
                
                OA = OscillationAnalyzer;
                
                CorResSU = [];
                cellCnt = 0;
                
                for iArea = 1:nAreas
                    
                    if isempty(results(iPatient).areasList{iArea})
                        continue
                    end
                    
                    % area = classifyArea(results(iPatient).areasList{iArea});
                    
                    % if isempty(area); error('need to update area list classification'); end
                    
                    %                     if whatTorun.singlePlots
                    %                         if exist('f0','var')
                    %                             set(f0, 'Position', get(0, 'Screensize'));
                    %                             saveas(f0,[folderToSaveFigures,'\',sprintf('%s.jpg',get(gcf,'name'))]);
                    %                             close(f0); clear f0
                    %                         end
                    %
                    %                         f0 = newA4figure(sprintf('cross_corr_%s_%s_1',results(iPatient).patientName,results(iPatient).areasList{iArea}));
                    %
                    %                         Mr = 5;
                    %                         fig_cnt = 1;
                    %                         fig_name_cnt = 1;
                    %
                    %                     end
                    
                    clear timestamps_vec
                    for iChannel = 1:length(channelsList{iArea})
                        for iCell = 1:size(results(iPatient).spikeCELL,3)
                            % Prepare spike vectors for cross-cor calculation
                            try
                                if      isempty(results(iPatient).spikeCELL{iArea,iChannel,iCell})
                                    continue
                                end
                            catch
                                continue
                            end
                            cellCnt = cellCnt + 1;
                            
                            f0 = newA4figure(sprintf('cross_corr_%s_ch%d_cell%d',results(iPatient).patientName,channelsList{iArea}(iChannel),iCell));
                            colorScheme = {'r','b','g','k'};
                            for ii_a = 1:4
                                COND = ii_a;
                                timestamps_vec{ii_a} = results(iPatient).spikeCELL{iArea,iChannel,iCell}.spikeTimes_ms{COND};
                                OA.TC_at_ROI = OA.CalcTemporalCrossCorrelation(timestamps_vec{ii_a},timestamps_vec{ii_a} );
                                CorResSU(cellCnt,ii_a,:) = OA.TC_at_ROI;
                                OA.power_intrinsic_TC_smoothed = OA.CalcPowerSpec(CorResSU(cellCnt,ii_a,:) );
                                power_intrinsic_TC_smoothed_SU(cellCnt,ii_a,:) = OA.power_intrinsic_TC_smoothed;
                                
                                if whatTorun.singlePlots
                                    % subplot(Mr,2,2*(fig_cnt-1)+1)
                                    subplot(3,1,1)
                                    hold on
                                    lineInfo.color = colorScheme{ii_a};
                                    cmap = [1 0 0; 0 1 0; 0 0 1; 0 0 0];
                                    OA.PlotCorrelogram(false, lineInfo)
                                    
                                    str = sprintf('auto-corr - single unit, area %s',results(iPatient).areasList{iArea});
                                    title(str,'fontsize',obj.crossCorfsz);
                                    axis tight
                                    XLIM = get(gca,'xlim');
                                    axis([XLIM 0 ,4*10^-3])
                                    
                                    
                                    if ii_a == 4
                                        axes('position',[0.7,0.75,0.15,0.15])
                                        samplingRate = 32e3;
                                        lsz = 8;
                                        fsz = 8;
                                        spikeInfo = results(iPatient).spikeCELL{iArea, iChannel,iCell}.spikeInfo;
                                        TT = [0:(1/samplingRate): (length(spikeInfo.spike_shapes_mean)-1)*(1/samplingRate)];
                                        hh = shadedErrorBar(TT ,...
                                            spikeInfo.spike_shapes_mean,spikeInfo.spike_shapes_std);
                                        YLIM = get(gca,'ylim');
                                        set(gca,'ylim',[YLIM(1),60])
                                        YLIM = get(gca,'ylim');
                                        XLIM = get(gca,'xlim');
                                        hold on
                                        plot([TT(end)-(1e-3),TT(end)],(YLIM(1)-.1*diff(YLIM))*ones(1,2),'linewidth',lsz,'color','k')
                                        text(TT(end)-0.5*diff(([TT(end)-(1e-3),TT(end)])),YLIM(1)-0.37*diff(YLIM),'1ms','horizontalalignment','center','fontsize',fsz)
                                        plot(TT(1)*ones(1,2)-0.1*diff(XLIM),[10 60],'linewidth',lsz,'color','k')
                                        text(TT(1)-diff(XLIM)*0.28,35,'50uV','horizontalalignment','center','fontsize',fsz,'rotation',90)
                                        axis([-1e-3 TT(end) (YLIM(1)-.1*diff(YLIM)) YLIM(2)])
                                        axis off;
                                    end
                                    
                                    
                                    subplot(3,1,2)
                                    hold on
                                    %subplot(Mr,2,2*(fig_cnt-1)+2)
                                    OA.PlotCorrelogramPowerSpec(false, lineInfo)
                                    legend(conditions)
                                    str = 'power spectrum of auto-corr';
                                    title(str,'fontsize',obj.crossCorfsz);
                                    
                                    % adding ERP reponses -
                                    subplot(3,1,3)
                                    hold all
                                    
                                    if ~strcmp(conditions{ii_a },'N')
                                        mm = matfile(fullfile('E:\MAZE\data_p\populationAnalysis\neuralUnits_trialData\',...
                                            sprintf('neuralUnitsTrials_%s_EXP%s_cond%s.mat',runData(iPatient).patientName,runData(iPatient).EXP,conditions{ii_a})));
                                        ratesStruct = mm.results;
                                        stimTriggeredRates = ratesStruct.stimTriggeredRates{iArea}{iChannel};
                                    else
                                        stimTriggeredRates = ratesStruct.contTriggeredRates{iArea}{iChannel}; % 'N is the control
                                    end
                                    
                                    
                                    TT = [-500:1:500];
                                    %                                     plot(TT,mean(stimTriggeredRates),'color',cmap(ii_a,:))
                                    
                                    hh = shadedErrorBar(TT ,...
                                        mean(stimTriggeredRates),std(stimTriggeredRates)/sqrt(size(stimTriggeredRates,1)));
                                    hh.mainLine.Color = cmap(ii_a,:);
                                    hh.mainLine.LineWidth = 2;
                                    axis tight
                                    YLIM = get(gca,'ylim');
                                    axis([-obj.corrBuffer obj.corrBuffer YLIM])
                                    xlabel('time (ms)' )
                                    ylabel('mean and sem (uV)')
                                    title('MU (all units on channel) activity')
                                end
                            end
                            
                            
                            
                            PrintActiveFigs(folderToSaveFigures)
                            
                        end % icell
                    end  % iChannel
                end % iArea
                
                
                % RESULTS
                results(iPatient).OscillationAnalyzerV = OA;
                results(iPatient).CorResSU = CorResSU;
                results(iPatient).power_intrinsic_TC_smoothed_SU = power_intrinsic_TC_smoothed_SU;
                
                
                if ~isempty(fileNameResults)
                    save(fileNameResults,'results','-v7.3');
                end
                
            end
            
            %             % MU results -
            %             for ii_a = 1:2
            %
            %                 if ii_a == 1
            %                     OA.TC_at_ROI = nanmean(crossCorRes_pre(logical(isFrontal_mu),:));
            %                     OA.power_intrinsic_TC_smoothed = OA.CalcPowerSpec(OA.TC_at_ROI);
            %                     lineInfo.color = 'b';
            %                 else
            %                     OA.TC_at_ROI = nanmean(crossCorRes_post(logical(isFrontal_mu),:));
            %                     OA.power_intrinsic_TC_smoothed = OA.CalcPowerSpec(OA.TC_at_ROI);
            %                     lineInfo.color = 'r';
            %                 end
            %
            %                 subplot(rn,cn,cn*(pop_fig_cnt-1)+3); hold on
            %                 OA.PlotCorrelogram(false, lineInfo)
            %                 title(sprintf('Mean cross-cor pt %s, frontal units',pt),'fontsize',obj.crossCorfsz)
            %                 axis tight
            %
            %
            %                 subplot(rn,cn,cn*(pop_fig_cnt-1)+4); hold on
            %                 OA.PlotCorrelogramPowerSpec(false, lineInfo)
            %                 title(sprintf('Mean power spec pt %s, frontal units',pt),'fontsize',obj.crossCorfsz)
            %                 legend('pre','post');
            %
            %             end
        end
        
        
        
        
        function MACROLFP_spectralAnalysis(obj, runData, outputFileFolder, outputFigureFolder, whatToRun)
            
            global data_p_path
            global BEHAV_TRIALS_FOLDER
            
            if nargin < 3
                fileNameResults = '';
            end
            
            nPatients = length(runData);
            
            
            %go over the patients requested
            for iPatient = 1:nPatients
                clear AREAS
                
                disp(['Patient ',runData(iPatient).patientName]);
                results.patientName = runData(iPatient).patientName;
                results.EXP = runData(iPatient).EXP;
                
                %load exp data for stimulation timings
                expData = load(runData(iPatient).ExpDataFileName);
                expData = expData.EXP_DATA;
                
                
                try
                    filename = dir(fullfile(BEHAV_TRIALS_FOLDER,sprintf('*%s*',results.patientName(2:end))));
                    
                    if strcmp(results.EXP, '2')
                        mm = matfile(fullfile(BEHAV_TRIALS_FOLDER, filename(2).name));
                    else
                        mm = matfile(fullfile(BEHAV_TRIALS_FOLDER, filename(1).name));
                    end
                    
                    eval(['data_structure = mm.da',results.patientName(4:end),'_data_structure;'])
                    test_data_structure = mm.test_data_structure;
                    
                    obj.samplingRate = data_structure.sampling_rate;
                    Nelectrodes = size(data_structure.D_rep1,3);
                catch
                    disp('lfp data not found')
                end
                
                try
                    disp(runData(iPatient).macroMontageFileName)
                    load(runData(iPatient).macroMontageFileName)
                catch
                    error('macro montage file not found')
                end
                
                for ii = 1:Nelectrodes
                    contact_id = data_structure.channels_included(ii);
                    AREAS{ii} = MacroMontage(contact_id).Area;
                end
                if strcmp(whatToRun.timeStampsType,'X')
                    idx2 = data_structure.X_rep2_good_maze_trials;
                    idx3 = data_structure.X_rep3_good_maze_trials;
                    trials = [data_structure.X_rep2(idx2,:,:); data_structure.X_rep3(idx3,:,:)];
                    idx2 = data_structure.N_rep2_good_maze_trials;
                    idx3 = data_structure.N_rep3_good_maze_trials;
                    cont_trials =  [data_structure.N_rep2(idx2,:,:); data_structure.N_rep3(idx3,:,:)];
                elseif strcmp(whatToRun.timeStampsType,'D')
                    idx2 = data_structure.D_rep2_good_maze_trials;
                    idx3 = data_structure.D_rep3_good_maze_trials;
                    trials = [data_structure.D_rep2(idx2,:,:); data_structure.D_rep3(idx3,:,:)];
                    idx2 = data_structure.N_rep2_good_maze_trials;
                    idx3 = data_structure.N_rep3_good_maze_trials;
                    cont_trials =  [data_structure.N_rep2(idx2,:,:); data_structure.N_rep3(idx3,:,:)];
                elseif strcmp(whatToRun.timeStampsType,'G')
                    idx2 = data_structure.G_rep2_good_maze_trials;
                    idx3 = data_structure.G_rep3_good_maze_trials;
                    trials = [data_structure.G_rep2(idx2,:,:); data_structure.G_rep3(idx3,:,:)];
                    idx2 = data_structure.N_rep2_good_maze_trials;
                    idx3 = data_structure.N_rep3_good_maze_trials;
                    cont_trials =  [data_structure.N_rep2(idx2,:,:); data_structure.N_rep3(idx3,:,:)];
                elseif strcmp(whatToRun.timeStampsType,'Gm')
                    idx2 = data_structure.G_rep2_good_maze_trials;
                    idx3 = data_structure.G_rep3_good_maze_trials;
                    trials = [data_structure.G_rep2(idx2,:,:); data_structure.G_rep3(idx3,:,:)];
                    cont_trials =   data_structure.G_rep1;
                elseif strcmp(whatToRun.timeStampsType,'Xm')
                    idx2 = data_structure.X_rep2_good_maze_trials;
                    idx3 = data_structure.X_rep3_good_maze_trials;
                    trials = [data_structure.X_rep2(idx2,:,:); data_structure.X_rep3(idx3,:,:)];
                    cont_trials =   data_structure.X_rep1;
                elseif strcmp(whatToRun.timeStampsType,'Test')
                    %DBG
                end
                
                
                windowLfpForComparisonBfrAftr_ms = floor(obj.windowLfpForComparison*obj.samplingRate/1000);
                midPoint = floor(size(trials,2)/2);
                indsForSign = [midPoint-windowLfpForComparisonBfrAftr_ms:midPoint+windowLfpForComparisonBfrAftr_ms];
                
                for iiC = 1:Nelectrodes
                    disp(sprintf('%d/%d electrodes',iiC,Nelectrodes))
                    %                     time_halfbandwidth = 2.5;
                    %                     dps_seq = dpss(size(trials,1),time_halfbandwidth);
                    %                     params.Fs = obj.samplingRate;
                    clear params
                    params.Fs = obj.samplingRate; % Hz
                    params.trialave  = true;
                    params.fpass = [0 150];
                    params.tapers=[3 5]; % very smoothed (default was [3 5]
                    p = 0.05;
                    params.err  = [2 p]; % Jackknife error bars
                    % estimate the spectrum of averaged single trials based on macro local field potential
                    % trials expected as times * trials
                    midpoint = size(trials,2)/2+1;
                    timerange = [midpoint-windowLfpForComparisonBfrAftr_ms: midpoint + windowLfpForComparisonBfrAftr_ms];
                    
                    [estSpectrum,freq,confidenceBands] = mtspectrumc( squeeze(trials(:,timerange,iiC))', params );
                    [estSpectrum_cont,freq,confidenceBands_cont] = mtspectrumc( squeeze(cont_trials(:,timerange,iiC))', params );
                    results.spectralAnalysis(iiC).estSpectrum = estSpectrum;
                    results.spectralAnalysis(iiC).estSpectrum_cont = estSpectrum_cont;
                    results.spectralAnalysis(iiC).freq = freq;
                    results.spectralAnalysis(iiC).confidenceBands = confidenceBands;
                    results.spectralAnalysis(iiC).confidenceBands_cont = confidenceBands_cont;
                    results.spectralAnalysis(iiC).spectrumCalcParams = params;
                    
                    % Add moving window to estimate the spectrum of averaged single trials based Ftaperson micro local field potential
                    % Note units have to be consistent. Thus, if movingwin is in seconds, Fs
                    % has to be in Hz. see chronux.m for more information.
                    clear params
                    params.Fs = obj.samplingRate; % Hz
                    params.fpass = [5 150];
                    params.trialave  = true;
                    params.tapers = [10 19];
                    movingwin = [0.5 0.05] ; % sec
                    params.pad = 0;
                    params.err=0;
                    [estSpectrum_movingwin,t_movingwin,freq_movingwin] = mtspecgramc( squeeze(trials(:,:,iiC))', movingwin, params );
                    [estSpectrum_movingwin_cont,t_movingwin,freq_movingwin] = mtspecgramc( squeeze(cont_trials(:,:,iiC))', movingwin, params );
                    
                    results.spectralAnalysis(iiC).estSpectrum_movingwin = estSpectrum_movingwin;
                    results.spectralAnalysis(iiC).estSpectrum_movingwin_cont = estSpectrum_movingwin_cont;
                    results.spectralAnalysis(iiC).t_movingwin = t_movingwin;
                    results.spectralAnalysis(iiC).freq_movingwin = freq_movingwin;
                    results.spectralAnalysis(iiC).spectrumCalcParams.params = params;
                    results.spectralAnalysis(iiC).spectrumCalcParams.movingwin = movingwin;
                    
                    
                    % Coherence between trials
                    clear params
                    params.Fs = obj.samplingRate; % Hz
                    params.fpass = [0 200];
                    params.trialave  = true;
                    params.tapers = [6 10];
                    params.pad = 0;
                    
                    % for iiC = 1:Nelectrodes
                    win = 0.1;
                    A = trials(:,:,iiC);
                    B = A';
                    try
                        [Sc,Cmat,Ctot,Cvec,Cent,f]=CrossSpecMatc(B,win,params);
                    catch
                        warning('Cmat Coherence failed')
                        Ctot = [];
                        f = [];
                    end
                    A = cont_trials(:,:,iiC);
                    B = A';
                    try
                        [Sc_cont,Cmat_cont,Ctot_cont,Cvec_cont,Cent_cont,f_cont]=CrossSpecMatc(B,win,params);
                    catch
                        warning('Cmat Coherence failed')
                        Ctot_cont = [];
                        f_cont = [];
                    end
                    results.spectralAnalysis(iiC).estTotalCoh = Ctot;
                    results.spectralAnalysis(iiC).estTotalCoh_cont = Ctot_cont;
                    results.spectralAnalysis(iiC).freq_coh = f;
                    % Output:
                    %       Sc (cross spectral matrix frequency x channels x channels)
                    %       Cmat Coherence matrix frequency x channels x channels
                    %       Ctot Total coherence: SV(1)^2/sum(SV^2) (frequency)
                    %       Cvec leading Eigenvector (frequency x channels)
                    %       Cent A different measure of total coherence: GM/AM of SV^2s
                    %       f (frequencies)
                    
                    
                    aucStim = sum(trials(:,indsForSign,iiC)');
                    aucCont = sum(cont_trials(:,indsForSign,iiC)');
                    try
                        [psA,~] = ranksum(aucStim,aucCont);
                    catch
                        pSA = 1;
                        warning('rank sum was not ran')
                    end
                    
                    results.LFP_stats(iiC).psA = psA;
                    
                    results.AREAS = AREAS;
                    results.channels_included = data_structure.channels_included;
                    results.trials = trials;
                    results.cont_trials = cont_trials;
                    
                    if whatToRun.newtimef
                        EEG.data1 = trials(:,:,iiC);
                        EEG.data2 = cont_trials(:,:,iiC);
                        EEG.pnts = length(EEG.data1);
                        MAXFREQ = 250;
                        
                        title_str_p = sprintf('elect_%d_%s',data_structure.channels_included(iiC),AREAS{iiC});
                        figureName = sprintf('%s_%s_MACRO_%s_newtimeF_%s',results.patientName,results.EXP, title_str_p, whatToRun.timeStampsType);
                        % Example to compare two condition (channel 1 EEG versus ALLEEG(2)):
                        %        >> [ersp,itc,powbase,times,freqs,erspboot,itcboot] = ...
                        %                  newtimef({EEG.data(1,:,:) ALLEEG(2).data(1,:,:)},
                        %                       EEG.pnts, [EEG.xmin EEG.xmax]*1000, EEG.srate, cycles);
                        % Outputs:
                        %            ersp   = (nfreqs,timesout) matrix of log spectral diffs from baseline
                        %                     (in dB log scale or absolute scale). Use the 'plot' output format
                        %                     above to output the ERSP as shown on the plot.
                        %            itc    = (nfreqs,timesout) matrix of complex inter-trial coherencies.
                        %                     itc is complex -- ITC magnitude is abs(itc); ITC phase in radians
                        %                     is angle(itc), or in deg phase(itc)*180/pi.
                        %          powbase  = baseline power spectrum. Note that even, when selecting the
                        %                     the 'trialbase' option, the average power spectrum is
                        %                     returned (not trial based). To obtain the baseline of
                        %                     each trial, recompute it manually using the tfdata
                        %                     output described below.
                        %            times  = vector of output times (spectral time window centers) (in ms).
                        %            freqs  = vector of frequency bin centers (in Hz).
                        %         erspboot  = (nfreqs,2) matrix of [lower upper] ERSP significance.
                        %          itcboot  = (nfreqs) matrix of [upper] abs(itc) threshold.
                        %           tfdata  = optional (nfreqs,timesout,trials) time/frequency decomposition
                        %                      of the single data trials. Values are complex.
                        
                        
                        [ersp,itc,powbase,times,freqs,erspboot,~] = ...
                            newtimef({EEG.data1', EEG.data2'},...
                            EEG.pnts , [1 1000], obj.samplingRate, 0, 'alpha',0.05, 'freqs' ,[0 MAXFREQ],'baseboot',0);
                        ff = gcf;
                        set(gcf,'name',figureName)
                        PrintActiveFigs(outputFigureFolder)
                        results.ersp(iiC,:,:) = ersp;
                        results.erspboot{iiC} = erspboot;
                        
                        figureName = sprintf('%s_%s_MACRO_w_baseline_%s_newtimeF_%s',results.patientName,results.EXP, title_str_p, whatToRun.timeStampsType);
                        [ersp,itc,powbase,times,freqs,erspboot,~] = ...
                            newtimef({EEG.data1', EEG.data2'},...
                            EEG.pnts , [-200 800], obj.samplingRate, 0, 'alpha',0.05, 'freqs' ,[0 MAXFREQ],'baseboot',0);
                        ff = gcf;
                        set(gcf,'name',figureName)
                        PrintActiveFigs(outputFigureFolder)
                    end
                end
                
                
                if ~isempty(outputFileFolder)
                    filename = sprintf('%s_E%s_stimuli_%s_MACRO_spectralAnalysis',runData(iPatient).patientName,results.EXP,whatToRun.timeStampsType);
                    fileNameResults = fullfile(outputFileFolder, filename);
                    save(fileNameResults,'results');
                end
                
                if whatToRun.singlePlots
                    if strcmp(whatToRun.timeStampsType,'GeGm')
                        ref_str = 'G exp';
                    else
                        ref_str = 'N';
                    end
                    
                    figureName = sprintf('%s_%s_MACRO_LFP_VALIDATION_%s',results.patientName,results.EXP, whatToRun.timeStampsType);
                    newA4figure(figureName);
                    r = ceil(sqrt(2*Nelectrodes+2));
                    for iiC = 1:Nelectrodes
                        subplot(r,r,2*(iiC-1)+1)
                        hold all
                        midPoint = size(trials,2)/2+1;
                        t = [1:size(trials(:,:,iiC),2)]/obj.samplingRate;
                        n = size(trials(:,:,iiC),1);
                        imagesc(t,1:n,trials(:,:,iiC))
                        axis ij
                        colorbar
                        
                        clear title_str
                        title_str{1} = sprintf('trials, elect %d, %s',data_structure.channels_included(iiC),AREAS{iiC});
                        title(title_str)
                        xlabel('t (ms)')
                        ylabel('trials')
                        
                        subplot(r,r,2*(iiC-1)+2)
                        hold all
                        n = size(cont_trials(:,:,iiC),1);
                        imagesc(t,1:n,cont_trials(:,:,iiC))
                        axis ij
                        colorbar
                        
                        title_str{1} = sprintf('control, elect %d, %s',data_structure.channels_included(iiC),AREAS{iiC});
                        title(title_str)
                        xlabel('t (ms)')
                        
                    end
                    
                    figureName = sprintf('%s_%s_MACRO_AverageLFP_%s',results.patientName,results.EXP, whatToRun.timeStampsType);
                    newA4figure(figureName);
                    M = 4;
                    if Nelectrodes > 16; M = 5; end
                    
                    for iiC = 1:Nelectrodes
                        subplot(M,M,iiC)
                        hold all
                        midpoint = size(trials,2)/2+1;
                        
                        Am = mean(trials(:,:,iiC));
                        As = std(trials(:,:,iiC))/sqrt(size(trials(:,:,iiC),2));
                        shadedErrorBar(t,Am,As, 'lineProps','r')
                        Am = mean(cont_trials(:,:,iiC));
                        As = std(cont_trials(:,:,iiC))/sqrt(size(cont_trials(:,:,iiC),1));
                        shadedErrorBar(t,Am,As)
                        idx2test = 1:length(Am);
                        starLoc = max(Am(:))*1.5;
                        clear signif
                        for iiT = 1:length(idx2test)
                            a = trials(:,idx2test(iiT),iiC);
                            b = cont_trials(:,idx2test(iiT),iiC);
                            
                            [T df] = ttest2_cell( { a b } ,'inhomogenous');
                            signif(iiT) = 2*tcdf(-abs(T), df(1));
                            if signif(iiT) < 0.05
                                plot(t(idx2test(iiT)),starLoc,'b*');
                            end
                        end
                        
                        clear title_str
                        title_str{1} = sprintf('electrode %d, %s',data_structure.channels_included(iiC),AREAS{iiC});
                        title_str{2} = sprintf('ttest %2.2e',results.LFP_stats(iiC).psA);
                        
                        plot(t(midPoint)*ones(1,2), get(gca,'ylim'),'k')
                        title(title_str)
                        xlabel('t (s)')
                        ylabel('uV')
                        if iiC == 1
                            legend(whatToRun.timeStampsType,ref_str)
                        end
                    end
                    
                    
                    figureName = sprintf('%s_%s_MACRO_spectralAnalysis1_%s',results.patientName,results.EXP, whatToRun.timeStampsType);
                    newA4figure(figureName);
                    for iiC = 1:Nelectrodes
                        subplot(M,M,iiC)
                        hold all
                        freq = results.spectralAnalysis(iiC).freq;
                        plot_vector(results.spectralAnalysis(iiC).estSpectrum, freq,1,...
                            results.spectralAnalysis(iiC).confidenceBands,'r')
                        plot_vector(results.spectralAnalysis(iiC).estSpectrum_cont, freq,1,...
                            results.spectralAnalysis(iiC).confidenceBands_cont,'k')
                        
                        ylabel('Spectrum dB')
                        
                        title_str = sprintf('electrode %d, %s, jackknife err, p =.05',data_structure.channels_included(iiC),AREAS{iiC});
                        title(title_str)
                        
                        P = get(gca,'position');
                        axes('Position',[P(1)+2*P(3)/4, P(2)+1.8*P(4)/4, P(3)/2, P(4)/2 ])
                        hold all
                        freq = results.spectralAnalysis(iiC).freq;
                        plot_vector(results.spectralAnalysis(iiC).estSpectrum, freq,1,...
                            results.spectralAnalysis(iiC).confidenceBands,'r')
                        plot_vector(results.spectralAnalysis(iiC).estSpectrum_cont, freq,1,...
                            results.spectralAnalysis(iiC).confidenceBands_cont,'k')
                        axis([0 50 0 inf])
                        ylabel('')
                        
                    end
                    
                    if(0) % this analysis is not working - need to work it out
                        figureName = sprintf('%s_%s_MACRO_spectralAnalysis2_%s',runData(iPatient).patientName, results.EXP,whatToRun.timeStampsType);
                        newA4figure(figureName);
                        for iiC = 1:Nelectrodes
                            subplot(M,M,iiC)
                            hold all
                            estSpectrum_movingwin =  mean(10*log10(abs(results.spectralAnalysis(iiC).estSpectrum_movingwin).^2),3);
                            cont_estSpectrum_movingwin =  mean(10*log10(abs(results.spectralAnalysis(iiC).estSpectrum_movingwin_cont).^2),3);
                            plot_diff =  estSpectrum_movingwin;
                            freq_movingwin = results.spectralAnalysis(iiC).freq_movingwin;
                            
                            a = plot_diff(:,5:end);
                            CLIM = [prctile(a(:),5) prctile(a(:),95) ];
                            h1 = imagesc(t_movingwin,freq_movingwin, plot_diff',CLIM);
                            axis xy
                            axis([t_movingwin(1) t_movingwin(end) 5 150])
                            set(gca,'xlim',[0.2 .8], 'XTick',[0.2 0.5 .8],'XTickLabel',{'-0.3', '0', '0.3'},'ytick',[5 50 100 150])
                            hold all
                            plot(0.5*ones(1,2),get(gca,'ylim'),'k')
                            colorbar
                            title(sprintf('trials - %s',whatToRun.timeStampsType))
                            ylabel('Frequency (Hz)')
                            xlabel('Time (sec)')
                            h = colorbar;
                            title_str = sprintf('electrode %d, %s',data_structure.channels_included(iiC),AREAS{iiC});
                            title(title_str)
                            
                            
                        end
                    end
                    
                    figureName = sprintf('%s_%s_MACRO_spectralAnalysis2a_%s',runData(iPatient).patientName, results.EXP,whatToRun.timeStampsType);
                    newA4figure(figureName);
                    suptitle('spectrum diff')
                    for iiC = 1:Nelectrodes
                        subplot(M,M,iiC)
                        hold all
                        estSpectrum_movingwin =  mean(10*log10(abs(results.spectralAnalysis(iiC).estSpectrum_movingwin).^2),3);
                        cont_estSpectrum_movingwin =  mean(10*log10(abs(results.spectralAnalysis(iiC).estSpectrum_movingwin_cont).^2),3);
                        plot_diff =  estSpectrum_movingwin-cont_estSpectrum_movingwin;
                        freq_movingwin = results.spectralAnalysis(iiC).freq_movingwin;
                        
                        a = plot_diff(:,5:end);
                        CLIM = [prctile(a(:),5) prctile(a(:),95) ];
                        h1 = imagesc(t_movingwin,freq_movingwin, plot_diff',CLIM);
                        axis xy
                        axis([t_movingwin(1) t_movingwin(end) 5 150])
                        set(gca,'xlim',[0.2 .8], 'XTick',[0.2 0.5 .8],'XTickLabel',{'-0.3', '0', '0.3'},'ytick',[5 50 100 150])
                        hold all
                        plot(0.5*ones(1,2),get(gca,'ylim'),'k')
                        colorbar
                        title(sprintf('trials - %s',whatToRun.timeStampsType))
                        ylabel('Frequency (Hz)')
                        xlabel('Time (sec)')
                        title_str = sprintf('electrode %d, %s',data_structure.channels_included(iiC),AREAS{iiC});
                        title(title_str)
                        h = colorbar;
                        h.Label.String= 'dB';
                        h.Label.Rotation =  270;
                        
                    end
                    
                    figureName = sprintf('%s_%s_MACRO_spectralAnalysis3_%s',runData(iPatient).patientName,results.EXP, whatToRun.timeStampsType);
                    newA4figure(figureName);
                    suptitle('total coherence per electrode')
                    for iiC = 1:Nelectrodes
                        subplot(M,M,iiC)
                        hold all
                        Ctot = results.spectralAnalysis(iiC).estTotalCoh;
                        if isempty(Ctot); continue; end
                        Ctot_cont = results.spectralAnalysis(iiC).estTotalCoh_cont;
                        if isempty(Ctot_cont); continue; end
                        f = results.spectralAnalysis(iiC).freq_coh;
                        title(sprintf('electrode %d',iiC));
                        hold all
                        plot(f, Ctot,'r'); hold all; plot(f, Ctot_cont,'k');
                        if iiC == 1
                            legend('trials',ref_str)
                        end
                        title_str = sprintf('electrode %d, %s',data_structure.channels_included(iiC),AREAS{iiC});
                        title(title_str)
                    end
                    
                    %                     CLIM = [prctile(10*log10(estSpectrum_movingwin(:)),15) prctile(10*log10(estSpectrum_movingwin(:)),85) ];
                    %                     h1 = imagesc(t_movingwin,freq_movingwin, 10*log10(estSpectrum_movingwin),CLIM);
                    %                     axis([0.1 0.5 0 200])
                    %                     set(gca,'xlim',[0.1 0.5], 'XTick',[0.1 0.3 0.5],'XTickLabel',{'-0.2', '0', '0.2'})
                    %                     hold all
                    %                     plot(0.3*ones(1,2),get(gca,'ylim'),'k')
                    %                     colorbar
                    %                     title(sprintf('trials - %s',whatToRun.timeStampsType))
                    %                     axis ij
                    %                     ylabel('Frequency (Hz)')
                    %                     xlabel('Time (sec)')
                    %
                    %                     subplot(3,2,4)
                    %                     h2 = imagesc(t_movingwin,freq_movingwin, 10*log10(estSpectrum_movingwin_cont),CLIM);
                    %                     title(sprintf('control trials - %s','N'))
                    %                     set(gca,'xlim',[0.1 0.5], 'XTick',[0.1 0.3 0.5],'XTickLabel',{'-0.2', '0', '0.2'})
                    %                     colorbar
                    %                     hold all
                    %                     plot(0.3*ones(1,2),get(gca,'ylim'),'k')
                    %                     xlabel('Time (sec)')
                    %                     axis ij
                    
                    PrintActiveFigs(outputFigureFolder)
                    
                end
                
            end % patients
            
            clear results
            
        end % function
        
        
        function MACROLFP_crossChannel_spectralAnalysis(obj, runData, outputFileFolder, outputFigureFolder, whatToRun)
            
            global data_p_path
            global BEHAV_TRIALS_FOLDER
            
            if nargin < 3
                fileNameResults = '';
            end
            
            nPatients = length(runData);
            
            
            %go over the patients requested
            for iPatient = 1:nPatients
                clear AREAS
                
                disp(['Patient ',runData(iPatient).patientName]);
                results.patientName = runData(iPatient).patientName;
                results.EXP = runData(iPatient).EXP;
                
                %load exp data for stimulation timings
                expData = load(runData(iPatient).ExpDataFileName);
                expData = expData.EXP_DATA;
                
                try
                    filename = dir(fullfile(BEHAV_TRIALS_FOLDER,sprintf('*%s*',results.patientName(2:end))));
                    
                    if strcmp(results.EXP, '2')
                        mm = matfile(fullfile(BEHAV_TRIALS_FOLDER, filename(2).name));
                    else
                        mm = matfile(fullfile(BEHAV_TRIALS_FOLDER, filename(1).name));
                    end
                    
                    eval(['data_structure = mm.da',results.patientName(4:end),'_data_structure;'])
                    obj.samplingRate = data_structure.sampling_rate;
                    Nelectrodes = size(data_structure.D_rep1,3);
                catch
                    disp('lfp data not found')
                end
                
                try
                    disp(runData(iPatient).macroMontageFileName)
                    load(runData(iPatient).macroMontageFileName)
                catch
                    error('macro montage file not found')
                end
                
                for ii = 1:Nelectrodes
                    contact_id = data_structure.channels_included(ii);
                    AREAS{ii} = MacroMontage(contact_id).Area;
                end
                if strcmp(whatToRun.timeStampsType,'X')
                    idx2 = data_structure.X_rep2_good_maze_trials;
                    idx3 = data_structure.X_rep3_good_maze_trials;
                    trials = [data_structure.X_rep2(idx2,:,:); data_structure.X_rep3(idx3,:,:)];
                    idx2 = data_structure.N_rep2_good_maze_trials;
                    idx3 = data_structure.N_rep3_good_maze_trials;
                    cont_trials =  [data_structure.N_rep2(idx2,:,:); data_structure.N_rep3(idx3,:,:)];
                elseif strcmp(whatToRun.timeStampsType,'D')
                    idx2 = data_structure.D_rep2_good_maze_trials;
                    idx3 = data_structure.D_rep3_good_maze_trials;
                    trials = [data_structure.D_rep2(idx2,:,:); data_structure.D_rep3(idx3,:,:)];
                    idx2 = data_structure.N_rep2_good_maze_trials;
                    idx3 = data_structure.N_rep3_good_maze_trials;
                    cont_trials =  [data_structure.N_rep2(idx2,:,:); data_structure.N_rep3(idx3,:,:)];
                elseif strcmp(whatToRun.timeStampsType,'G')
                    idx2 = data_structure.G_rep2_good_maze_trials;
                    idx3 = data_structure.G_rep3_good_maze_trials;
                    trials = [data_structure.G_rep2(idx2,:,:); data_structure.G_rep3(idx3,:,:)];
                    idx2 = data_structure.N_rep2_good_maze_trials;
                    idx3 = data_structure.N_rep3_good_maze_trials;
                    cont_trials =  [data_structure.N_rep2(idx2,:,:); data_structure.N_rep3(idx3,:,:)];
                elseif strcmp(whatToRun.timeStampsType,'Gm')
                    idx2 = data_structure.G_rep2_good_maze_trials;
                    idx3 = data_structure.G_rep3_good_maze_trials;
                    trials = [data_structure.G_rep2(idx2,:,:); data_structure.G_rep3(idx3,:,:)];
                    cont_trials =   data_structure.G_rep1;
                end
                
                windowLfpForComparisonBfrAftr_ms = floor(obj.windowLfpForComparison*obj.samplingRate/1000);
                midPoint = floor(size(trials,2)/2);
                indsForSign = [midPoint-windowLfpForComparisonBfrAftr_ms:midPoint+windowLfpForComparisonBfrAftr_ms];
                
                
                if whatToRun.isMTLPt
                    
                    if isempty(runData(iPatient).HPC_channels_spectralAnalysis) || ...
                            isempty(runData(iPatient).OFC_channels_spectralAnalysis); continue; end
                    
                    HipCh = runData(iPatient).HPC_channels_spectralAnalysis(1);
                    FrontalCh = runData(iPatient).OFC_channels_spectralAnalysis(1);
                    
                    
                    ind1 = find(ismember(data_structure.channels_included, HipCh));
                    ind2 = find(ismember(data_structure.channels_included, FrontalCh));
                    A  = trials(:,:,ind1);
                    A_cont  = cont_trials(:,:,ind1);
                    B  = trials(:,:,ind2);
                    B_cont  = cont_trials(:,:,ind2);
                    EEG.icaact(1,:) = A(:);
                    ALLEEG.icaact(1,:) = A_cont(:);
                    EEG.icaact(2,:) = B(:);
                    ALLEEG.icaact(2,:) = B_cont(:);
                    
                    EEG.pnts = 1000;
                    MAXFREQ = 250;
                    
                    title_str_p = sprintf('elect_%d_%s_%d_%s',data_structure.channels_included(ind1),AREAS{ind1},data_structure.channels_included(ind2),AREAS{ind2});
                    figureName = sprintf('%s_%s_MACRO_%s_newcrossfF_%s',results.patientName,results.EXP, title_str_p, whatToRun.timeStampsType);
                    
                    
                    [coh,mcoh,timesout,freqsout,cohboot,cohangles,...
                        allcoher,alltfX,alltfY] = newcrossf({EEG.icaact(1,:) ...
                        ALLEEG.icaact(1,:)},{EEG.icaact(2,:) ...
                        ALLEEG.icaact(2,:)},EEG.pnts,[-200 800],obj.samplingRate, ...
                        0, 'alpha',0.05, 'freqs' ,[0 MAXFREQ], 'alpha',nan);
                    % results  =
                    ff = gcf;
                    set(gcf,'name',figureName)
                    PrintActiveFigs(outputFigureFolder)
                end
                
                
                if ~isempty(outputFileFolder)
                    filename = sprintf('%s_E%s_stimuli_%s_MACRO_crossCh_spectralAnalysis',runData(iPatient).patientName,results.EXP,whatToRun.timeStampsType);
                    fileNameResults = fullfile(outputFileFolder, filename);
                    save(fileNameResults,'results');
                end
            end
            
        end
        
        function AMNA_population_analysis(obj, runData, outputFileFolder, outputFigureFolder, whatToRun)
            
            global data_p_path;
            global MACROLFP_spectralAnalysisFolder;
            pVal_TH = 0.05;
            
            if nargin < 3
                fileNameResults = '';
            end
            nPatients = length(runData);
            
            cnt = 1;
            electrodeID_ALL = [];
            clear ptId  AREA_ALL
            %go over the patients requested
            for iPatient = 1:nPatients
                clear pVal_AUC
                disp(['Patient ',runData(iPatient).patientName]);
                results_all(iPatient).patientName = runData(iPatient).patientName;
                results_all(iPatient).EXP = runData(iPatient).EXP;
                
                try
                    filename = sprintf('%s_E%s_stimuli_%s_MACRO_spectralAnalysis',runData(iPatient).patientName,...
                        runData(iPatient).EXP,whatToRun.timeStampsType);
                    
                    fileNameResults = fullfile(outputFileFolder, filename);
                    
                    mm = matfile(fileNameResults)';
                    results = mm.results;
                    
                    macroMontage = load(runData(iPatient).macroMontageFileName);
                catch
                    error('data file not found')
                end
                
                try
                    load(runData(iPatient).macroMontageFileName)
                catch
                    error('macromontage not found')
                end
                
                nelectrodes = size(results.channels_included,2);
                
                % calculate coherence based on bands
                bands(1,:) = [1 4];
                bands(2,:) = [4 8];
                bands(3,:) = [8 12];
                bands(4,:) = [15 30];
                bands(5,:) = [30 200];
                
                clear coherence_sum coherenceCont_sum pVal_AUC
                for iiC = 1:nelectrodes
                    for ii_b = 1:size(bands)
                        
                        Ctot = results.spectralAnalysis(iiC).estTotalCoh;
                        Ctot_cont = results.spectralAnalysis(iiC).estTotalCoh_cont;
                        
                        f = results.spectralAnalysis(iiC).freq_coh;
                        if ~isempty(f) && ~isempty(Ctot_cont)
                            range(1) = find(f<=bands(ii_b,1),1,'last');
                            if isempty(find(f>=bands(ii_b,2),1,'first'))
                                range(2) = length(f);
                            else
                                range(2) = find(f>=bands(ii_b,2),1,'first');
                            end
                            
                            coherence_sum(iiC,ii_b) = sum(Ctot(range,:));
                            coherenceCont_sum(iiC,ii_b) = sum(Ctot_cont(range,:));
                        else
                            coherence_sum(iiC,ii_b) = nan;
                            coherenceCont_sum(iiC,ii_b) = nan;
                        end
                        pVal_AUC(iiC) = results.LFP_stats(iiC).psA;
                    end
                    
                    ptId{cnt} = results.patientName;
                    electrodeID_ALL = [electrodeID_ALL, results.channels_included(iiC)];
                    AREA_ALL{cnt}  = results.AREAS{iiC};
                    diffCoherence_ALL{cnt} = (coherence_sum - coherenceCont_sum)./(coherence_sum + coherenceCont_sum);
                    
                    LFP_average_DIFF(cnt) = results.LFP_stats(iiC).psA < pVal_TH;
                    
                    cnt = cnt+1;
                end
                
                results_all(iPatient).coherence_sum = coherence_sum;
                results_all(iPatient).coherenceCont_sum = coherenceCont_sum;
                results_all(iPatient).diffCoherence = (coherence_sum - coherenceCont_sum)./(coherence_sum + coherenceCont_sum);
                results_all(iPatient).psA = pVal_AUC;
                results_all(iPatient).AREAS = results.AREAS;
                results_all(iPatient).channels = results.channels_included;
                
                clear results
                
            end % patients
            
            if ~isempty(outputFileFolder)
                filename = sprintf('populationSummaru_spectralAnalysis_%s',whatToRun.timeStampsType);
                fileNameResults = fullfile(outputFileFolder, filename);
                save(fileNameResults, 'results_all','ptId','electrodeID_ALL','AREA_ALL','diffCoherence_ALL','LFP_average_DIFF')
            end
        end %pop analysis
        
        
    end % methods
    
end % class


