addpath('C:\Users\mgeva\Documents\GitHub\external\KK_useful\useful\violin')
CMAP = brewermap(3,'Dark2');

newA4figure('summary_AUC_Ratios');
set(gcf,'DefaultAxesFontSize',16);

markerSize = 16;
cnt = 0;
fileDir = 'E:\MAZE\data_p\populationAnalysis\macroLFP_ERP\';

for plotLine = 1:3

    if plotLine == 1
        rep = '2';
        DATAx = TABLE_ALL.aucCI2;
        DATAy = TABLE_ALL.aucCI1;
        ylabel_plot = 'Ratio rep 2/1';
        xlabel_plot = 'MTL/Frontal contacts';
        subplotID = 1:3;

    elseif plotLine == 2
        rep = '3';
        DATAx = TABLE_ALL.aucCI3;
        DATAy = TABLE_ALL.aucCI1;
        ylabel_plot = 'Ratio rep 3/1';
        subplotID = 4:6;

    elseif plotLine == 3
        DATAx = TABLE_ALL.aucCIMem;
        DATAy = TABLE_ALL.aucCI1;
        ylabel_plot = 'Ratio Memory/Exporation';
        subplotID = 7:9;


    elseif plotLine == 4
        DATAx = TABLE_ALL.aucContrastAcrossRep1Rep2;
        DATAy = TABLE_ALL.aucContrastAcrossRep1Rep3;
        ylabel_plot = 'auc_1 - auc_3/auc_1 + auc_3';
    end

    for ii_a = 1:3


        if ii_a == 1
            load(fullfile(fileDir,'allPts_stimuli_X_MACRO_aucTable.mat'))
            str = 'Planning';
        elseif ii_a == 2
            load(fullfile(fileDir,'allPts_stimuli_G_MACRO_aucTable.mat'))
            str = 'Goal';
        else
            load(fullfile(fileDir,'allPts_stimuli_D_MACRO_aucTable.mat'))
            str = 'Decision';
        end
     
        if plotLine == 1
            rep = '2';
            DATAx = TABLE_ALL.aucCI2;
            DATAy = TABLE_ALL.aucCI1;
            ylabel_plot = 'Ratio rep 2/1';
            subplotID = 1:2;

        elseif plotLine == 2
            rep = '3';
            DATAx = TABLE_ALL.aucCI3;
            DATAy = TABLE_ALL.aucCI1;
            ylabel_plot = 'Ratio rep 3/1';
            subplotID = 4:5;

        elseif plotLine == 3
            DATAx = TABLE_ALL.aucCIMem;
            DATAy = TABLE_ALL.aucCI1;
            ylabel_plot = 'Ratio Memory/Exporation';
            subplotID = 7:8;

        end

        subplot(3,3,subplotID) 
        hold all;

        % [XX,TFrm1] = rmoutliers(DATAx);
        % [YY,TFrm2] = rmoutliers(DATAy);
        % legalIDs = ~ismember( [1:length(DATAx)],unique([find(TFrm1(:)'),find(TFrm2(:)')]));
        % DATAx = DATAx(legalIDs);
        % DATAy = DATAy(legalIDs);
        RatioDD = DATAx./DATAy;
        [XX,TFrm1] = rmoutliers(RatioDD);
        legalIDs = ~ismember( [1:length(RatioDD)],find(TFrm1(:)));
        RatioDD = RatioDD(legalIDs);

        loc = ii_a;
        linecolor = CMAP(ii_a,:);
        lw = [10 1];
        my_boxplot(RatioDD,loc,linecolor,lw)
        [p(ii_a),h] = signrank(RatioDD);
        axis([0.5 3.5 -6 6])
        hold all;

        plot(get(gca,'xlim'),[0 0],'k')
        index = find(TABLE_ALL.isMTL(legalIDs));
        scatterX = loc - 0.2 + rand([1,length(index)])*0.4;
        hold all; plot(scatterX, RatioDD(index), 'm.','MarkerSize',markerSize)
        N(1) = length(index);

        index = find(TABLE_ALL.isHip(legalIDs));
        scatterX =loc - 0.2 + rand([1,length(index)])*0.4;
        hold all; plot(scatterX, RatioDD(index), 'r.','MarkerSize',markerSize)
        N(2) = length(index);

        index = find(TABLE_ALL.isFrontal(legalIDs));
        scatterX = loc - 0.2+ rand([1,length(index)])*0.4;
        hold all; plot(scatterX, RatioDD(index), 'b.','MarkerSize',markerSize)
        N(3) = length(index);

        xlabel_str{ii_a} = sprintf('%s_{p} = %2.2e,%d,%d,%d',str,p(ii_a),N(1),N(2),N(3))

        if ii_a == 3
            set(gca,'xtick',1:3,'XTickLabel',xlabel_str,'fontsize',7, 'ytick',[-6 0 6])
            title(sprintf('Rep %s/Rep1:%s',rep))
            ylabel(ylabel_plot)
        end

    end
end
saveas(gcf, fullfile('E:\MAZE\data_p\populationAnalysis\macroLFP_spectralAnalysis\figures\',[get(gcf,'name'),'.tif']))

%% seperate population from specific areas
markerSize = 12;
cnt = 0;

for ii_a = 1:3
    if ii_a == 1
        load(fullfile(fileDir,'allPts_stimuli_X_MACRO_aucTable.mat'))
        str = 'Planning';

    elseif ii_a == 2
        load(fullfile(fileDir,'allPts_stimuli_G_MACRO_aucTable.mat'))
        str = 'Goal';
    else
        load(fullfile(fileDir,'allPts_stimuli_D_MACRO_aucTable.mat'))
        str = 'Decision';
    end
    newA4figure(sprintf('%s',str))

    for plotLine = 1:3

        if plotLine == 1
            rep = '2';
            DATAx = TABLE_ALL.aucCI2;
            DATAy = TABLE_ALL.aucCI1;
            plotStr = sprintf('auc rep 1,2 CI (%s vs. Neutral)',str);

        elseif plotLine == 2
            rep = '3';
            DATAx = TABLE_ALL.aucCIMem;
            DATAy = TABLE_ALL.aucCI1;
            plotStr = sprintf('auc mem vs. 1 CI (%s vs. Neutral)',str);

        elseif plotLine == 3
            DATAx = TABLE_ALL.aucContrastAcrossRep1Rep2;
            DATAy = TABLE_ALL.aucContrastAcrossRep1Rep3;
            plotStr = sprintf('%s - Rep1 vs. Rep2, Rep1 vs. Rep3',str);
        end

        subplot(3,4,(plotLine-1)*4 + 1)
        histogram(DATAx,[-1.5:0.1:1.5])
        hold all
        histogram(DATAy,[-1.5:0.1:1.5])
        title(sprintf('%s',plotStr))

        cmap = brewermap(8,'dark2');
        subplot(3,4,(plotLine-1)*4 + 2)
        index = find(TABLE_ALL.isMTL);
        clear A
        A(1,:) = DATAy(index);
        A(2,:) = DATAx(index);
        [p(1),h] = signrank(A(1,:),A(2,:));

        %violinplot(1,A(1,:)','color',cmap(1,:),'cutoff',0.1)
        my_boxplot(A(1,:)',1,cmap(1,:))
        hold all
        % violinplot(2,A(2,:)','color',cmap(2,:))
        my_boxplot(A(2,:)',2,cmap(2,:))

        clear A
        index = find(TABLE_ALL.isFrontal);
        A(1,:) = DATAy(index);
        A(2,:) = DATAx(index);
        [p(2),h] = signrank(A(1,:),A(2,:));
        %         violinplot(3,A(1,:)','color',cmap(1,:))
        %         hold all
        %         violinplot(4,A(2,:)','color',cmap(2,:))
        my_boxplot(A(1,:)',3,cmap(1,:))
        hold all
        % violinplot(2,A(2,:)','color',cmap(2,:))
        my_boxplot(A(2,:)',4,cmap(2,:))
        set(gca,'xtick',[1:4],'xticklabels',{'mtl-rep1','mtl-rep2','frontal-rep1','frontal-rep2'},'XTickLabelRotation',45)
        title(sprintf('p_m = %2.2e,p_f = %2.2e',p(1),p(2)))

        if plotLine == 1 & ii_a == 2
            set(gca,'ytick',[-4:4:4])
            axis([0.5 4.5 -4 4])
        end
        plot(get(gca,'xlim'),zeros(1,2),'k')


        subplot(3,4,(plotLine-1)*4 + 3)
        hold all
        index = find(TABLE_ALL.isMTL);
        clear A
        A = DATAx(index);
        [p(1),h] = signrank(A);
        my_boxplot(A(:)',1,cmap(1,:))

        clear A
        index = find(TABLE_ALL.isFrontal);
        A = DATAx(index);
        [p(2),h] = signrank(A);
        my_boxplot(A(:)',2,cmap(2,:))

        set(gca,'xtick',[1:4],'xticklabels',{'MTLCh','FrontalCh'})
        axis([0.5 2.5 -4 4])

        title(sprintf('MTLp = %2.2e, Fp = %2.2e',p(1),p(2)))
        plot(get(gca,'xlim'),zeros(1,2),'k')


    end
end

%%
figure
markerSize = 12;
cnt = 0;
for plotLine = 1:3
    for ii_a = 1:3
        if ii_a == 1
            load(fullfile(fileDir,'allPts_stimuli_X_MACRO_aucTable.mat'))
            str = 'Planning';
        elseif ii_a == 2
            load(fullfile(fileDir,'allPts_stimuli_G_MACRO_aucTable.mat'))
            str = 'Goal';
        else
            load(fullfile(fileDir,'allPts_stimuli_D_MACRO_aucTable.mat'))
            str = 'Decision';
        end
        cnt = cnt + 1;
        if plotLine == 1
            rep = '2';

            DATAx = TABLE_ALL.aucRatio2;
            DATAy = TABLE_ALL.aucRatio1;
            ylabel_plot = 'ratio rep 1';
            xlabel_plot = 'ratio auc rep 2';
        elseif plotLine == 2
            rep = '3';

            DATAx = TABLE_ALL.aucRatio3;
            DATAy = TABLE_ALL.aucRatio1;
            ylabel_plot = 'ratio rep 1';
            xlabel_plot = 'ratio auc rep 3';
        elseif plotLine == 3
            rep = '1 vs rep';

            DATAx = TABLE_ALL.aucRatioAcrossRep1Rep2;
            DATAy = TABLE_ALL.aucRatioAcrossRep1Rep3;
            ylabel_plot = 'auc ratio rep 1/rep2';
            xlabel_plot = 'auc ratio rep1/rep 3';
        end
        subplot(3,3,cnt)
        plot(DATAx, DATAy, 'k.','MarkerSize',markerSize)
        axis([-300 300 -300 300])
        hold all;
        [p(ii_a),h] = signrank(DATAx, DATAy);

        %axis([-1 1 -1 1])
        %axis([-1 1 -1 1]*2)
        %hold all; plot([-1 1]*2,[-1 1]*2,'--k')

        plot(get(gca,'xlim'),[0 0],'k')
        plot([0 0],get(gca,'ylim'),'k')
        set(gca,'xtick',[-300 0 300],'ytick',[-300 0 300])
        title(sprintf('Rep1 vs Rep %s:%s, %2.2e',rep,str,p(ii_a)))


        index = find(TABLE_ALL.isMTL);
        hold all; plot(DATAx(index), DATAy(index), 'm.','MarkerSize',markerSize)

        index = find(TABLE_ALL.isHip);
        hold all; plot(DATAx(index), DATAy(index), 'r.','MarkerSize',markerSize)

        index = find(TABLE_ALL.isFrontal);
        hold all; plot(DATAx(index),DATAy(index), 'b.','MarkerSize',markerSize)

        xlabel('REP2')
        ylabel('REP1')
    end
end


%% Collecting data to a table for modeling work


%%
p_improve = [];
for bb = 1:2

    bType = {'G','X'};
    fileList = dir(fullfile(fileDir, sprintf('*_behavNeural_%s*',bType{bb})));
    MTLChan_auc_all = [];
    RT_diff_rep1_rep2_MTL_all = [];
    frontalChan_auc_all  = [];
    RT_diff_rep1_rep2_frontal_all  = [];
    cmap = brewermap(30,'Dark2');
    for ii_a = 1
        mtl_cnt = 1;
        f_cnt= 1;

        if ii_a == 1
            filename = sprintf('TABLE_compare_rep1AUC_toBehav_%s',bType{bb});
        else
            filename = sprintf('TABLE_compare_rep1AUC_toBehav3_%s',bType{bb});
        end

        for iP = 1:length(fileList)

            mm = matfile(fullfile(fileDir, fileList(iP).name));
            p_improve(bb,iP) = mm.p_improve;

            if sum(mm.P_corr < 0.05)
                disp(fileList(iP).name)
                disp('specific ch correlated with behav')
            end

            mazeCompletionTable = mm.mazeCompletionTable;
            missingMazes = mm.missingMazes;

            if ii_a == 1
                NeuralMarker = mm.aucStimR1_Ch_mazes;
                BehavMarker = mm.RT_diff_rep1_rep2;
            else
                NeuralMarker = mm.aucStimR1_Ch_mazes;
                BehavMarker = mm.RT_diff_rep2_rep3;
            end

            MTLChan_auc = (NeuralMarker(logical(mm.MTLCh),:));
            MTLChan_auc(abs(MTLChan_auc) == 1) = NaN;
            frontalChan_auc = (NeuralMarker(logical(mm.frontCh),:));
            frontalChan_auc(abs(frontalChan_auc) == 1) = NaN;

            MTL_ids = find(mm.MTLCh);
            if ~isempty(MTL_ids)
                for ii = 1:sum(mm.MTLCh)

                    ptID = ones(1,length(BehavMarker))*iP;
                    NeuralD = MTLChan_auc(ii,:);
                    ChID = ones(1,length(BehavMarker))*mm.channelsIncluded(MTL_ids(ii),1);
                    isM = ones(1,length(BehavMarker));
                    isF = zeros(1,length(BehavMarker));
                    TablePerCh = table(ptID(:), NeuralD(:),ChID(:), isM(:),isF(:),BehavMarker(:),...
                        'VariableNames',{'pt','neuralMarker','chID','isMTL','isFrontal','behavMarker'});

                    if mtl_cnt == 1 & f_cnt == 1
                        TABLE_ALL = TablePerCh;
                    else

                        TABLE_ALL = [TABLE_ALL; TablePerCh];
                    end
                    mtl_cnt= mtl_cnt + 1;

                end
            end
            F_ids = find(mm.frontCh);
            if ~isempty(F_ids)
                for ii = 1:sum(mm.frontCh)

                    ptID = ones(1,length(BehavMarker))*iP;
                    NeuralD = frontalChan_auc(ii,:);
                    ChID = ones(1,length(BehavMarker))*mm.channelsIncluded(F_ids(ii),1);
                    isF = ones(1,length(BehavMarker));
                    isM = zeros(1,length(BehavMarker));
                    TablePerCh = table(ptID(:), NeuralD(:),ChID(:), isM(:),isF(:),BehavMarker(:),...
                        'VariableNames',{'pt','neuralMarker','chID','isMTL','isFrontal','behavMarker'});

                    if mtl_cnt == 1 & f_cnt == 1
                        TABLE_ALL = TablePerCh;
                    else

                        TABLE_ALL = [TABLE_ALL; TablePerCh];
                    end
                    f_cnt= f_cnt + 1;

                end
            end

            if ~isempty(frontalChan_auc)
                try
                    mdl = fitlm(frontalChan_auc',BehavMarker(:),'RobustOpts','on');
                    ss = anova(mdl,'summary');
                    f_ch_p(iP) = ss.pValue(2);
                catch
                    disp('not enough data')
                    continue
                end

                %         set(0,'CurrentFigure',fig1)
                %         subplot(3,3,iP)
                %         plot(mdl)
                %         title(sprintf('%s - ses %d, P = %2.2e',bType{ii}, iP,  f_ch_p(iP)))
                %         suptitle('frontal channels')
            end

            if ~isempty(MTLChan_auc)
                try
                    mdl = fitlm(MTLChan_auc',BehavMarker(:),'RobustOpts','on');
                    ss = anova(mdl,'summary');
                    m_ch_p(iP) = ss.pValue(2);
                catch
                    disp('not enough data')
                    continue
                end
            end
            clear frontalChan_auc MTLChan_auc
        end
        behavMarkerZ = zscore(TABLE_ALL.behavMarker);
        TABLE_ALL.behavMarkerZ = behavMarkerZ;

        outputFileFolder = fullfile(data_p_path,'populationAnalysis\macroLFP_ERP');
        save(fullfile(outputFileFolder,filename),'TABLE_ALL','runData')
        clear TABLE_ALL
    end
end


%%
%Zscore the behavioral marker - which is the time-diff between the
%completion of a maze on the first try and the second try
outputFileFolder = fullfile(data_p_path,'populationAnalysis\macroLFP_ERP');

for bb = 1:2

    bType = {'G','X'};
    disp(bType{bb})

    filename = sprintf('TABLE_compare_rep1_rep2_%s',bType{bb});
    mm = matfile(fullfile(outputFileFolder,filename));
    TABLE_ALL = mm.TABLE_ALL;

    idx1 = find(TABLE_ALL.isMTL);
    formula = 'behavMarkerZ ~ 1 + (1|pt)+(1|chID)';
    lme1 = fitlme(TABLE_ALL(idx1,:),formula);
    formula = 'behavMarkerZ ~ neuralMarker + (1|pt)+(1|chID)';
    lme2 = fitlme(TABLE_ALL(idx1,:),formula);
    results = compare(lme1,lme2);
    disp('mtl neural marker predicts behav')
    disp(results.pValue)

    TABLE_ALL_M = TABLE_ALL(idx1,:);
    pts = unique(TABLE_ALL_M.pt);
    for ii = 1:length(pts)
        ses = pts(ii);
        idx1 = (TABLE_ALL_M.pt == pts(ii));
        MM{ii} = (TABLE_ALL_M.neuralMarker(idx1));
    end

    idx1 = find(TABLE_ALL.isFrontal);
    formula = 'behavMarkerZ ~ 1 + (1|pt)+(1|chID)';
    lme1 = fitlme(TABLE_ALL(idx1,:),formula);
    formula = 'behavMarkerZ ~ neuralMarker + (1|pt)+(1|chID)';
    lme2 = fitlme(TABLE_ALL(idx1,:),formula);
    results = compare(lme1,lme2);
    disp('frontal neural marker predicts behav')
    disp(results.pValue)

end