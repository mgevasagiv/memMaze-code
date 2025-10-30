outputFileFolder = fullfile(data_p_path,'populationAnalysis\macroLFP_ERP');

for bb = 1:2    
    bType = {'G','X'};    
    for ii_a = 1:3
        mtl_cnt = 1;
        f_cnt= 1;
        
        if ii_a == 1
            filename = sprintf('TABLE_compare_rep1_rep2_%s.mat',bType{bb});
            str = 'rep1 vs rep 2';
        elseif ii_a == 2
            filename = sprintf('TABLE_compare_rep1_rep3_%s.mat',bType{bb});
            str = 'rep2 vs rep 3';
        elseif ii_a == 3
            filename = sprintf('TABLE_compare_rep1AUC_toBehav_%s.mat',bType{bb});
            str = 'AUC on 1';
        end
        load(fullfile(outputFileFolder,filename))
        disp(sprintf('linear mixed models testing MTL and frontal ch - %s %s\n',bType{bb},str))
        
        
        %%
        idx1 = find(TABLE_ALL.isMTL);
        formula = 'behavMarkerZ ~ 1 + (1|pt)+(1|chID)';
        lme1 = fitlme(TABLE_ALL(idx1,:),formula);
        formula = 'behavMarkerZ ~ neuralMarker + (1|pt)+(1|chID)';
        lme2 = fitlme(TABLE_ALL(idx1,:),formula);
        results = compare(lme1,lme2);
        disp(sprintf('MTL - P = %2.4e\n',results.pValue(2)))

        idx1 = find(TABLE_ALL.isFrontal);
        formula = 'behavMarkerZ ~ 1 + (1|pt)+(1|chID)';
        lme1 = fitlme(TABLE_ALL(idx1,:),formula);
        formula = 'behavMarkerZ ~ neuralMarker + (1|pt)+(1|chID)';
        lme2 = fitlme(TABLE_ALL(idx1,:),formula);
        results = compare(lme1,lme2);
        disp(sprintf('Frontal - P = %2.2f\n',results.pValue(2)))  
    end
end


%% ------------------
newA4figure('allPts')
idx1 =  find(TABLE_ALL.isMTL);
TABLE = TABLE_ALL(idx1,:);
% rmv outliers
idx3_rmv = abs(TABLE.behavMarker) > 3*std(abs(TABLE.behavMarker));
TABLE(idx3_rmv,:) = [];

idx1 = TABLE.behavMarkerZ > 0;
idx2 = TABLE.behavMarkerZ < 0;
hold all
cla
my_boxplot(zscore(TABLE.neuralMarker(idx1)),1,'k')
my_boxplot(zscore(TABLE.neuralMarker(idx2)),2,'k')

formula = 'behavMarkerZ ~ neuralMarker + (1|pt)';
lme = fitlme(TABLE,formula,'FitMethod','REML');
beta = fixedEffects(lme);
[~,~,STATS] = randomEffects(lme);

TABLE = TABLE_ALL;
idx3_rmv = abs(TABLE.behavMarker) > 100;
TABLE(idx3_rmv,:) = [];
lm = fitlm(TABLE,'linear','RobustOpts','on');

