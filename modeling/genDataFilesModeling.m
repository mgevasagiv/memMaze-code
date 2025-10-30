% Generating files which allow matching brain activity to behav
codePath_MAZE();

fileDir = 'E:\MAZE\data_p\populationAnalysis\macroLFP_ERP\';
bType = {'G','X'};
cmap = brewermap(30,'Dark2');

p_improve = [];
for bb = 1:2
    
    fileList = dir(fullfile(fileDir, sprintf('*_behavNeural_%s*',bType{bb})));
    MTLChan_auc_all = [];
    RT_diff_rep1_rep2_MTL_all = [];
    frontalChan_auc_all  = [];
    RT_diff_rep1_rep2_frontal_all  = [];
    for ii_a = 1:2
        mtl_cnt = 1;
        f_cnt= 1;
        
        if ii_a == 1        
            filename = sprintf('TABLE_compare_rep1_rep2_%s',bType{bb});
        else
            filename = sprintf('TABLE_compare_rep1_rep3_%s',bType{bb});
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
                NeuralMarker = mm.aucContrastAcrossRep1Rep2_allMazes;
                BehavMarker = (mm.RT_diff_rep1_rep2);
                % mazeIDs = mazeCompletionTable.maze_id;
            else
                NeuralMarker = mm.aucContrastAcrossRep2Rep3_allMazes;
                BehavMarker = mm.RT_diff_rep2_rep3;
                % mazeIDs = mazeCompletionTable.maze_id;
            end
            id_rep = mm.id_rep;
            mazeIDs = id_rep{1};

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
                    TablePerCh = table(ptID(:), NeuralD(:),ChID(:), isM(:),isF(:),mazeIDs(:), BehavMarker(:),...
                        'VariableNames',{'pt','neuralMarker','chID','isMTL','isFrontal','mazeIDs','behavMarker'});
                    
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
                    
                    if length(mazeIDs) ~= length(NeuralD)
                        disp('1')
                    end
                    TablePerCh = table(ptID(:), NeuralD(:),ChID(:), isM(:),isF(:),mazeIDs(:),BehavMarker(:),...
                        'VariableNames',{'pt','neuralMarker','chID','isMTL','isFrontal','mazeIDs','behavMarker'});
                    
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
        pTListId = table([1:length(patients)]',patients',expNames','VariableNames',{'ptID','patients','expNames'});
        save(fullfile(outputFileFolder,filename),'TABLE_ALL','runData','pTListId')
        clear TABLE_ALL
    end
end