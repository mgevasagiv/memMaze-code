bands(1,:) = [1 4];
bands(2,:) = [4 8];
bands(3,:) = [8 12];
bands(4,:) = [15 30];
bands(5,:) = [30 200];



newA4figure('populationCoherenceChange')
for ii_C = 1:4
    clear A
    
    if ii_C == 1
        load('E:\MAZE\data_p\populationAnalysis\macroLFP_spectralAnalysis\populationSummaru_spectralAnalysis_X.mat')
        titlestr = 'First step (planning), vs Neutral points';
    elseif ii_C == 2
        load('E:\MAZE\data_p\populationAnalysis\macroLFP_spectralAnalysis\populationSummaru_spectralAnalysis_D.mat')
        titlestr = 'Decision making (memory based)';
    elseif ii_C == 3
        load('E:\MAZE\data_p\populationAnalysis\macroLFP_spectralAnalysis\populationSummaru_spectralAnalysis_G.mat')
        titlestr = 'Goal (memory based vs neutral)';
    elseif ii_C == 4
        load('E:\MAZE\data_p\populationAnalysis\macroLFP_spectralAnalysis\populationSummaru_spectralAnalysis_Gm.mat')
        titlestr = 'Goal (memory based vs exploration)';
    end
    
    subplot(2,2,ii_C)
    A =cell(1,5);
    AREA_STR = [];
    CMAP = brewermap(length(results_all),'Dark2');
    for iiS = 1:length(results_all)
        color_ses = CMAP(iiS,:);
        diffCoherence = results_all(iiS).diffCoherence;
        %         if iiS == 2
        %             diffCoherence = diffCoherence(end,:);
        %         end
        Nses = length(diffCoherence);
        
        cnt = length(AREA_STR);
        for iis_a = 1:Nses
            cnt = cnt+1;
            AREA_STR{cnt} =  results_all(iiS).AREAS{iis_a};
        end
        
        clear vecB
        auc_sig = results_all(iiS).psA(1:Nses) <0.05;
        
        for iib = 1:size(diffCoherence,2)
            vecB = diffCoherence(:,iib);
            if isnan(vecB); continue; end
            hold all
            scatterV = randn(1,Nses)*0.1;
            plot(iib*ones(1,Nses)+scatterV, vecB,'.','color',color_ses);
            plot(iib*ones(1,sum(auc_sig))+scatterV(auc_sig), vecB(auc_sig),'*','color',color_ses);
            A{iib} = [A{iib} vecB(:)'];
        end
        if Nses ~= length(results_all(iiS).AREAS)
            warning('%s missing contacts', results_all(iiS).patientName)
        end
    end
    
    text_st = {'pt1','pt2','pt3','pt3','pt4','pt5','pt6'};
    
    for iib = 1:5
        [p(ii_C,iib), hh(ii_C, iib)] = signrank(A{iib});
    end
    set(gca,'xtick',[1:5],'xticklabels',{'[0,4]','[4,8]','[8 12]','[15 30]','[30 200]'},'xlim',[0 6]...
        ,'ylim',[-0.45 0.45])
    plot(get(gca,'xlim'),[0 0],'k')
    title(titlestr)
    
    
    
    % Check whether classification can be performed based on a specific
    % frequecy
    coherence_sum_agg = [];
    coherenceCont_sum_agg = [];
    for ii = 1:length(results_all)
        
        coherence_sum_agg = [coherence_sum_agg;  results_all(ii).coherence_sum];
        coherenceCont_sum_agg = [coherenceCont_sum_agg;  results_all(ii).coherenceCont_sum];
    end
    
    %% Look at single-trial decoding based on different channels
    for iif = 1:5
        Nchannels = size(coherence_sum_agg,1);
        D1 = [coherence_sum_agg(:,iif); coherenceCont_sum_agg(:,iif)];
        gm = fitgmdist(D1,2,'RegularizationValue',0.000005);
        P = posterior(gm,D1);
        
        C1 = find(P(:,1)>P(:,2));
        C2 = find(P(:,2)>P(:,1));
        Svec = zeros(1,length(D1));
        if gm.mu(1) > gm.mu(2)
            Svec(C1) = 1;
        else
            Svec(C2) = 1;
        end
        Phit(iif) = sum(Svec(1:Nchannels))/Nchannels;
        Pfa(iif) = sum(Svec(Nchannels+1:end))/Nchannels;
    end
    
    net = feedforwardnet(100);
    D1 = [coherence_sum_agg; coherenceCont_sum_agg]';
    targets = [ones(1,Nchannels),zeros(1,Nchannels)];
    net = configure(net,D1,targets);
    [net,tr] = train(net,D1,targets);
    y1 = net(D1);
    [z b] = sort((sum(abs(net.IW{1}))));
    
    plot(y1,'.')
    hold all
    
    plot(get(gca,'xlim'),0.5*ones(1,2))
    Phit_NN = nansum(y1(1:Nchannels))/Nchannels;
    Pfa_NN = nansum(y1(Nchannels+1:end))/Nchannels;
    
    str{1} = sprintf('Phit = %1.2f ',Phit_NN);
    str{2} = sprintf('Pfa = %1.2f',Pfa_NN);
    str{3} = sprintf('highest weight - %d (%2.0f%%) %d (%2.0f%%)',b(5),z(b(5))*100/sum(z),b(4),z(b(4))*100/sum(z));
    XLIM = get(gca,'xlim');
    YLIM = get(gca,'ylim');
    text(XLIM(1),YLIM(1)-diff(YLIM)/4.5,str)
    
    
    
end

%% Looking at freq bands
load('E:\MAZE\data_p\populationAnalysis\macroLFP_spectralAnalysis\populationSummaru_spectralAnalysis_X.mat')
freq_movingwin = results_all(1).freq_movingwin;

clear p
for iiF = 1:6
    if iiF == 1
        evalSTR = 'Gamma_power_change_PRE_dB';
        newA4figure('populationGammaBandPreChange');
        freq_plot_id = find(freq_movingwin > 40);
    elseif iiF == 2
        evalSTR = 'Gamma_power_change_POST_dB';
        newA4figure('populationGammaBandPostChange');
        freq_plot_id = find(freq_movingwin > 40);
    elseif iiF == 3
        evalSTR = 'Ripple_power_change_PRE_dB';
        newA4figure('Ripple_power_change_PRE_dB');
        freq_plot_id = find(freq_movingwin > 40);
    elseif iiF == 4
        evalSTR = 'Ripple_power_change_POST_dB';
        newA4figure('Ripple_power_change_POST_dB');
        freq_plot_id = find(freq_movingwin > 40);
    elseif iiF == 5
        evalSTR = 'th_power_change_PRE_dB';
        newA4figure('th_power_change_PRE_dB');
        freq_plot_id = find(freq_movingwin < 40);
    elseif iiF == 6
        evalSTR = 'th_power_change_POST_dB';
        newA4figure('th_power_change_POST_dB');
        freq_plot_id = find(freq_movingwin < 40);
    end
    
    
    for ii_C = 1:4
        clear A
        diffGamma_all_F = [];
        diffGamma_all_M = [];
        if ii_C == 1
            load('E:\MAZE\data_p\populationAnalysis\macroLFP_spectralAnalysis\populationSummaru_spectralAnalysis_X.mat')
            titlestr = 'First step (planning), vs Neutral points';           
        elseif ii_C == 2
            load('E:\MAZE\data_p\populationAnalysis\macroLFP_spectralAnalysis\populationSummaru_spectralAnalysis_D.mat')
            titlestr = 'Decision making (memory based)';
        elseif ii_C == 3
            load('E:\MAZE\data_p\populationAnalysis\macroLFP_spectralAnalysis\populationSummaru_spectralAnalysis_G.mat')
            titlestr = 'Goal (memory based vs neutral)';
        elseif ii_C == 4
            load('E:\MAZE\data_p\populationAnalysis\macroLFP_spectralAnalysis\populationSummaru_spectralAnalysis_Gm.mat')
            titlestr = 'Goal (memory based vs exploration)';
        end
        cntM = 1;
        cntF = 1;
        spec_all_M = [];
        spec_all_F = [];
        
        subplot(4,4,(ii_C-1)*4+1)
        hold all
        A =cell(1,5);
        AREA_STR = [];
        CMAP = brewermap(length(results_all),'Dark2');
        for iiS = 1:length(results_all)
            color_ses = CMAP(iiS,:);
                            
            diffGammaZscore = (eval(['results_all(iiS).',evalSTR]));

            MTLCh = logical(results_all(iiS).MTLCh);
            if sum(MTLCh)
                diffGamma = diffGammaZscore(MTLCh);
                scatterV = randn(1,length(diffGamma))*0.1;
                plot(1*ones(1,length(diffGamma))+scatterV,zscore(diffGamma),'.','color',color_ses);
                diffGamma_all_M = [diffGamma_all_M, diffGamma];
                spec_all_M(cntM:cntM+sum(MTLCh)-1,:,:) = results_all(iiS).estSpectrum_movingwin_diff_dB(MTLCh,:,:);
                cntM = cntM+sum(MTLCh);   
            end
            
            fCh = logical(results_all(iiS).frontCh);
            if sum(fCh)
                diffGamma = diffGammaZscore(fCh);
                scatterV = randn(1,length(diffGamma))*0.1;
                plot(2*ones(1,length(diffGamma))+scatterV,diffGamma,'.','color',color_ses);
                diffGamma_all_F = [diffGamma_all_F, diffGamma];
                spec_all_F(cntF:cntF+sum(fCh)-1,:,:) = results_all(iiS).estSpectrum_movingwin_diff_dB(fCh,:,:);
                cntF = cntF+sum(fCh);
                
            end
            
            set(gca,'XTick',[1 2], 'XTickLabel',{'Hip/Am','prefr ctx'},'XLim',[0.5 2.5])
            ylabel('dB')
            
        end
        %     text_st = {'pt1','pt2','pt3','pt3','pt4','pt5','pt6'};
        %
        %     for iib = 1:5
        %         [p(ii_C,iib), hh(ii_C, iib)] = signrank(A{iib});
        %     end
        %     set(gca,'xtick',[1:5],'xticklabels',{'[0,4]','[4,8]','[8 12]','[15 30]','[30 200]'},'xlim',[0 6]...
        %         ,'ylim',[-0.45 0.45])
        [p(ii_C,1),h] = signrank(diffGamma_all_M);
        [p(ii_C,2),h] = signrank(diffGamma_all_F);
     
        
        plot(get(gca,'xlim'),[0 0],'k')
        title_str_all{1} = titlestr;
        title_str_all{2} = sprintf('p_{MTL} = %2.2f, p_{F} = %2.2f',p(ii_C,1),p(ii_C,2));
        title(title_str_all)
        
           
        subplot(4,4,(ii_C-1)*4+2)
        hold all
        bar([1 2],[mean((diffGamma_all_M)),mean((diffGamma_all_F))],'k')
        errorbar([1 2],[mean((diffGamma_all_M)),mean((diffGamma_all_F))],[std((diffGamma_all_M))/sqrt(length(diffGamma_all_M)),...
                            std((diffGamma_all_F)/sqrt(length(diffGamma_all_F)))],'k.')
        set(gca,'XTick',[1 2], 'XTickLabel',{'Hip/Am','prefr ctx'},'XLim',[0.5 2.5])

        
        freq_movingwin = results_all(iiS).freq_movingwin;
        t_movingwin = results_all(iiS).t_movingwin;
        subplot(4,4,(ii_C-1)*4+3)
        imagesc(t_movingwin, freq_movingwin(freq_plot_id), squeeze(mean(spec_all_M(:,:,(freq_plot_id))))')
        axis xy
        colorbar
        title('MTL spectral average)')
        hold all
        plot(ones(1,2)*0.5,get(gca,'ylim'),'k')
        
        subplot(4,4,(ii_C-1)*4+4)
        imagesc(t_movingwin, freq_movingwin(freq_plot_id), squeeze(mean(spec_all_F(:,:,(freq_plot_id))))')
        colorbar
        hold all
        plot(ones(1,2)*0.5,get(gca,'ylim'),'k')
        axis xy
        title('Frontal spectral average)')
        
        diffFreqBand{ii_C}.diffGamma_all_M = diffGamma_all_M;
        diffFreqBand{ii_C}.diffGamma_all_F = diffGamma_all_F;
        diffFreqBand{ii_C}.spec_all_M = spec_all_M;
        diffFreqBand{ii_C}.spec_all_F = spec_all_F;
        
    end
    
end

%% 
outputFigureFolder = 'E:\Dropbox\RanganathLab\MAZE\POSTER_IMAGES';
newA4figure('spectralAverage')
set(gcf,'DefaultAxesFontSize',24);

ii_C = 3; % Goal based on memory 
axes('position',[0.1 0.1 0.25 0.25])
imagesc(t_movingwin, freq_movingwin(19:end), squeeze(mean(diffFreqBand{ii_C}.spec_all_M (:,:,(19:end))))')

hold all
plot(ones(1,2)*0.5,get(gca,'ylim'),'k')
set(gca,'xtick',[0.25 0.5 0.75],'ytick',[60 100 140])
axis xy
title('MTL spectral profile - Goal vs N')
xlabel('t(sec)')
ylabel('f (hz)')
c = colorbar;
c.Ticks = [-1.2 0 0.6];

ii_C = 3; % Goal based on memory 
axes('position',[0.55 0.1 0.25 0.25])
imagesc(t_movingwin, freq_movingwin(19:end), squeeze(mean(diffFreqBand{ii_C}.spec_all_F (:,:,(19:end))))')

set(gca,'xtick',[0.25 0.5 0.75],'ytick',[60 100 140])
axis xy
xlabel('t(sec)')
ylabel('f (hz)')
title('Frontal spectral profile - Goal vs N')

c = colorbar;
c.Ticks = [-0.25 0 0.75];
hold all

plot(ones(1,2)*0.5,get(gca,'ylim'),'k')


ii_C = 4; % 
axes('position',[0.1 0.5 0.25 0.25])
imagesc(t_movingwin, freq_movingwin(19:end), -squeeze(mean(diffFreqBand{ii_C}.spec_all_M (:,:,(19:end))))')

hold all
plot(ones(1,2)*0.5,get(gca,'ylim'),'k')
set(gca,'xtick',[0.25 0.5 0.75],'ytick',[60 100 140])
axis xy
title('MTL spectral profile - 1st Goal vs Mem')
xlabel('t(sec)')
ylabel('f (hz)')
c = colorbar;
c.Ticks = [-1.2 0 0.6];

ii_C = 4; % 
axes('position',[0.55 0.5 0.25 0.25])
imagesc(t_movingwin, freq_movingwin(19:end), -squeeze(mean(diffFreqBand{ii_C}.spec_all_F (:,:,(19:end))))')

set(gca,'xtick',[0.25 0.5 0.75],'ytick',[60 100 140])
axis xy
xlabel('t(sec)')
ylabel('f (hz)')
title('Frontal spectral profile - 1st Goal vs Mem')

c = colorbar;
c.Ticks = [-0.25 0 0.75];
hold all

plot(ones(1,2)*0.5,get(gca,'ylim'),'k')


PrintActiveFigs(outputFigureFolder)