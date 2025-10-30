% Creating example - single channel - prediction of behavior (hip channel)
outputFileFolder = fullfile(data_p_path,'populationAnalysis\macroLFP_ERP');
fileDir = outputFileFolder;
filename = sprintf('TABLE_compare_rep1_rep2_G');
load(fullfile(outputFileFolder,filename))

pt = 3;
newA4figure('Example - mtl responses')
set(gcf,'DefaultAxesFontSize',24);
ax = axes('position',[0.1 0.3 0.12 0.2]);
idx1 =  find(TABLE_ALL.isMTL & (TABLE_ALL.pt == pt));
TABLE = TABLE_ALL(idx1,:);
idx1 = TABLE.behavMarkerZ > 0;
idx2 = TABLE.behavMarkerZ < 0;
hold all
[p,h] = ranksum(TABLE.neuralMarker(idx1), TABLE.neuralMarker(idx2));
my_boxplot(TABLE.neuralMarker(idx1),1,'k')
my_boxplot(TABLE.neuralMarker(idx2),2,'k')
set(gca,'XTick',[1 2],'XTickLabel',{'Improved','Degraded'},'XTickLabelRotation',35,'ytick',[-4 0 4])
axis([0.5 2.5 -4 4])
plot(get(gca,'xlim'),[0 0],'k')
plot([1 2],[4.2 4.2],'k','linewidth',2)
text(0.5, -10, sprintf('p = %2.2e',p))

% outputFigureFolder = 'E:\Dropbox\RanganathLab\MAZE\POSTER_IMAGES';
% PrintActiveFigs(outputFigureFolder)
%%
%% TBD - switching to box plot and maybe line plots for hippocampal channels
outputFigureFolder = fullfile(data_p_path,'populationAnalysis\macroLFP_ERP\figures');
sx = 0.12;
sy = 0.2;

% left to right, top to bottom
pos(1,:) = [0.1 0.75 sx sy];
pos(2,:) = [0.3 0.75 sx sy];
pos(3,:) = [0.1 0.5 sx sy];
pos(4,:) = [0.3 0.5 sx sy];
pos(5,:) = [0.1 0.25 sx sy];
pos(6,:) = [0.3 0.25 sx sy];
for iiF = 1:2
    if iiF == 1
        newA4figure(sprintf('summary_auc_poster',str))
        set(gcf,'DefaultAxesFontSize',24);
        noTitle = true;
    else
        newA4figure(sprintf('summary_auc_poster_stats',str))
        set(gcf,'DefaultAxesFontSize',10);
        noTitle = false;
    end
    for ii_a =1:6
        
        
        if ii_a == 1
            
            load(fullfile(fileDir,'allPts_stimuli_G_MACRO_aucTable.mat'))
            str = 'Goal';
            DATAx = TABLE_ALL.aucCI1;
            tstr{1} = 'auc goal, first rep';
        elseif ii_a == 2
            
            load(fullfile(fileDir,'allPts_stimuli_G_MACRO_aucTable.mat'))
            str = 'Goal';
            DATAx = TABLE_ALL.aucCIMem;
            tstr{1} = 'auc goal, suc. Mem mazes';
        elseif ii_a == 3
            
            load(fullfile(fileDir,'allPts_stimuli_X_MACRO_aucTable.mat'))
            str = 'X';
            DATAx = TABLE_ALL.aucCI1;
            tstr{1} = 'auc X, first rep';
        elseif ii_a == 4
            
            load(fullfile(fileDir,'allPts_stimuli_X_MACRO_aucTable.mat'))
            str = 'X';
            DATAx = TABLE_ALL.aucCIMem;
            tstr{1} = 'auc X, suc. Mem mazes';
            
            
        elseif ii_a == 5
            
            load(fullfile(fileDir,'allPts_stimuli_D_MACRO_aucTable.mat'))
            str = 'Dec';
            DATAx = TABLE_ALL.aucCI1;
            tstr{1} = 'auc dec, first rep';
        elseif ii_a == 6
            
            load(fullfile(fileDir,'allPts_stimuli_D_MACRO_aucTable.mat'))
            str = 'Dec';
            DATAx = TABLE_ALL.aucCIMem;
            tstr{1} = 'auc dec, suc. Mem mazes';
        end
        
        axes('position',pos(ii_a,:))
        
        hold all
        index = find(TABLE_ALL.isMTL);
        clear A
        A = DATAx(index);
        B = A;
        [p(1),h] = signrank(A);
        my_boxplot(A(:)',1,'k')
        
        axis([0.5 2.5 -4 5])
        ax = get(gca);
        ax.Clipping = "off";          
        if p(1) < 0.05
            plot([0.8 1.2],[4.5,4.5],'k','linewidth',3)
            % plot(1,5,'k*','markersize',6)
        end
        clear A
        index = find(TABLE_ALL.isFrontal);
        A = DATAx(index);
        [p(2),h] = signrank(A);
        my_boxplot(A(:)',2,[0.3 0.3 .9])
                ax = get(gca);
        ax.Clipping = "off";          
        if p(2) < 0.05
            plot([1.8 2.2],[4.5,4.5],'k','linewidth',3)
            % plot(1,5,'k*','markersize',6)
        end
        if sum(ii_a == [5,6])
            set(gca,'xtick',[1:2],'xticklabels',{'MTL','Frontal'},'XTickLabelRotation',35)
        else
            set(gca,'xtick',[5:6])
        end
        axis square
        set(gca,'ytick',[-4 0 4])
        tstr{2} = sprintf('MTLp = %2.2e, Fp = %2.2e',p(1),p(2));
        if ~noTitle
            title(tstr)
        end
        plot(get(gca,'xlim'),zeros(1,2),'k')

        [pp(ii_a) h] = ranksum(A(:),B(:));
        
        
        clear p
    end
end
        
% PrintActiveFigs(outputFigureFolder)

%% Adding behavior
fileDir  = 'E:\MAZE\data_p\populationAnalysis';
behav_table = readtable(fullfile(fileDir, 'behaviorSummary.xlsx'));
outputFigureFolder = 'E:\MAZE\data_p\populationAnalysis\behavFigures';

newA4figure('behavSummary')
set(gcf,'DefaultAxesFontSize',24);
xS = 0.2;
yS = 0.2;
CMAP = brewermap(7,'Dark2');
for iiT = 2:3
    
    if iiT == 1
        idx = find(strcmp(behav_table.RT_type, 'rtX'));
        plotTable = behav_table(idx,:);
        nP = size(plotTable,1)/3;
        axes('position',[0.1,0.15,xS,yS])
    elseif iiT == 2
        idx = find(strcmp(behav_table.RT_type, 'rtN'));
        plotTable = behav_table(idx,:); 
        axes('position',[0.35,0.15,xS,yS])
    elseif iiT == 3
        idx = find(strcmp(behav_table.RT_type, 'rtD'));
        plotTable = behav_table(idx,:);
        axes('position',[0.65,0.15,xS,yS])
    end
    
    
    
    for iiP = 1:7
        hold all
        ptD = plotTable.RT((iiP-1)*3 + 1 : (iiP-1)*3 + 3);
        ptD(2) = ptD(2)/ptD(1);
        ptD(3) = ptD(3)/ptD(1);
        ptD(1) = 1;
        plot(1:3,ptD,'.-','color',CMAP(iiP,:),'markersize',14,'linewidth',3);
        
        allP_RT(iiT,iiP,:) =  plotTable.RT((iiP-1)*3 + 1 : (iiP-1)*3 + 3);
        clear ptD
    end
    axis([0.5 6.5 0.7 1.05])
    set(gca,'xtick',[1:3],'XTickLabel',{'1','2','3'},'XTickLabelRotation',45)
    set(gca,'ytick',[0.7,1])
    
    if iiT == 2
        ylabel('RT (Norm)')
    end
    
    PP = get(gca,'position');
    axes('position',[PP(1)+PP(3)*0.85,PP(2)+0.65*PP(3),xS/2,yS*0.75])
    hold all
    my_boxplot(allP_RT(iiT,:,1),1,'k')
    my_boxplot(allP_RT(iiT,:,3),2,'k')
    set(gca,'xtick',[1:2],'XTickLabel',{'1','3'},'XTickLabelRotation',35)
    %if iiT == 2
        ylabel('RT (mS)')
    %end
    
    if sum(iiT == [3,2])
        axis([0.5 2.5 0 500])
        ax = gca;
        ax.Clipping = "off";
        set(gca,'ytick',[0, 500])
        plot(1:2,[500,500],'k','linewidth',3)
                plot(1.5,570,'k*','markersize',6)

        
    else
        axis([0.5 2.5 0 2000])
        set(gca,'ytick',[0, 2000])
        
    end
    
end

for iiT = 2:3
    [p(iiT),h(iiT)] = ttest(allP_RT(iiT,:,1),allP_RT(iiT,:,3));
end

% PrintActiveFigs(outputFigureFolder)
    %%
outputFileFolder = fullfile(data_p_path,'populationAnalysis\macroLFP_ERP');
data_p_path = 'E:\MAZE\Data_p\';
filename = sprintf('TABLE_compare_rep1_rep2_G');
mm = matfile(fullfile(outputFileFolder,filename));
disp(sprintf('linear mixed models testing MTL and frontal ch - %s %s\n',bType{bb},str))
       
TABLE_ALL = mm.TABLE_ALL;
pt = unique(TABLE_ALL.pt);
idx3_rmv = abs(TABLE_ALL.behavMarker) > 100;
TABLE_ALL(idx3_rmv,:) = [];

newA4figure('NeuralBehav_hist')
set(gcf,'DefaultAxesFontSize',24);

hold all
CMAP = brewermap(length(pt),'Dark2');
clf
for ii_a = 1:2
    if ii_a == 1
        axes('position',[0.1 0.2 0.15 0.15])
        histC = 'k';
    else
        axes('position',[0.4 0.2 0.15 0.15])
        histC = 'b';
    end
    hold all
    cnt = 0;
    clear M1 M2 M1B M2B
    for iP = 1:length(pt);
        if ii_a == 1
            idx1 =  find(TABLE_ALL.isMTL & (TABLE_ALL.pt == pt(iP)));
        else
            idx1 =  find(TABLE_ALL.isFrontal & (TABLE_ALL.pt == pt(iP)));
        end
        if isempty(idx1); continue; end
        TABLE = TABLE_ALL(idx1,:);
        
        ch = unique(TABLE.chID);
        
        for iiC = 1:length(ch)
            cnt = cnt + 1;
            idx1 =  (TABLE.pt == pt(iP)) & (TABLE.chID == ch(iiC)) ;
            TABLE_CH = TABLE(idx1,:);
            
            idx1 = TABLE_CH.behavMarker > 0;
            idx2 = TABLE_CH.behavMarker < 0;
            M1(cnt) = nanmean(TABLE_CH.neuralMarker(idx1));
            if isnan(M1(cnt))
                disp('1')
            end
            M2(cnt) = nanmean(TABLE_CH.neuralMarker(idx2));
            if isnan(M2(cnt))
                disp('1')
            end
            M1B(cnt) = nanmean(TABLE_CH.behavMarkerZ(idx1));
            M2B(cnt) = nanmean(TABLE_CH.behavMarkerZ(idx2));
            
%             plot([M1(cnt) M2(cnt)],[M1B(cnt) M2B(cnt)],'-','color',CMAP(iP,:))
%             plot([M1(cnt) ], [M1B(cnt) ],'o','color','k')
%             plot([M2(cnt) ], [M2B(cnt) ],'X','color','k')
        end

    end
%     axis([-1.5,2.5,  -1.2,0.5])
%     xlabel(sprintf('Neural change \n(Norm)'))
%     ylabel('\DeltaTime to finish (Zscore)')
%     set(gca,'ytick',[-1 0 0.5])
%     plot(get(gca,'xlim'),[0 0],'k','LineWidth',1)
%     
%     POS = get(gca,'position');
%     axes('position',[POS(1)+0.65*POS(3) POS(2)+0.8*POS(4) POS(3)*2/3 POS(4)*2/3])
    hold all
    vec = M1-M2;
    histogram(vec,[-1:0.25:1], 'Normalization','probability','FaceColor',histC,'EdgeColor',histC)
    axis([-1.5,1.5,  0 0.3])
    plot([0 0],get(gca,'ylim'),'k','LineWidth',1.5)
    set(gca,'ytick',[0 0.3])
    ylabel('Prob')
    % xlabel(sprintf('Neural change \n(Norm)'))
    
    [p h] = signtest(vec)
    VEC{ii_a} = vec;
    plot([1 1]*median(vec),get(gca,'ylim'),'r','LineWidth',2.5)
end

outputFigureFolder = 'E:\Dropbox\RanganathLab\MAZE\POSTER_IMAGES';
PrintActiveFigs(outputFigureFolder)
