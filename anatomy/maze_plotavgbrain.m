global globalFsDir;
globalFsDir='E:\MAZE\anatomy\FreesurferSubjectFolders';

patients = {'da017','da018', 'da019','da019','da022','da023','ir103'};
patients_id_anat = {'sub-017','sub-018','sub-019','sub-019_2','sub-022','sub-023','ir103'};
behav_str_vec = {'G','Gm','D','X','Xm','Test'};
summarYfilesFolder = 'E:\MAZE\data_p\populationAnalysis\macroLFP_spectralAnalysis\';

%%
PgroupAvgCoords=[];
PgroupLabels=[];
PgroupIsLeft=[];
NEURAL_MARKERS = [];
NEURAL_MARKERS_COHERENCE = [];

ptID = [];
ptList = {};

% Load markers to plot on averaged brain
for ii_a = 1:6
    
    behav_str = behav_str_vec{ii_a};
    fullfile(summarYfilesFolder, sprintf('populationSummaru_spectralAnalysis_%s.mat',behav_str));
    mm = matfile(summaryFile);
    results_all = mm.results_all;
    
    missingContacts = [];
    probe_elec_coord = cell(1,1);
    contactID = [];
    
    for iiP = 1:6
        currPt = results_all(iiP).patientName;
        disp(currPt)
        ptList{iiP} = currPt;
        results_pt = results_all(iiP);
        
        % Load electrode coordinates in native space
        anat_pt = patients_id_anat{iiP};
        disp(anat_pt)
        
        subPath = fullfile(globalFsDir,char(anat_pt));
        elecReconPath=fullfile(subPath,'elec_recon');
        filename = fullfile(elecReconPath, sprintf('%sPostimpLoc.txt',char(anat_pt)));
        [elec_name, elec_n, x, y, z, Hem, D] = textread(filename,'%s %d %f %f %f %s %s', 300);
        
        
        % Get the full contact name for each contact
        
        % load macro montage for area name
        macroMontageFileName = fullfile('E:\MAZE\data_p\MACRO_MONTAGE\', sprintf('%s',currPt),'EXP1','MacroMontage.mat');
        if strcmp(anat_pt,'sub-019_2')
            macroMontageFileName = fullfile('E:\MAZE\data_p\MACRO_MONTAGE\', sprintf('%s',currPt),'EXP2','MacroMontage.mat');
        end
        
        macroMontage = load(macroMontageFileName);
        macroMontage = macroMontage.MacroMontage;
        
        for iContact = 1:length(results_pt.channels)
            ChNum = results_pt.channels(iContact);
            ChArea = macroMontage(ChNum).Area;
            cnt = macroMontage(ChNum).Channel;
            ChLabel = sprintf('%s%d',ChArea,cnt);
            if ChNum == 1
                ChLabel = sprintf('%s%d',ChArea,ChNum);
            else
                c = true; cnt = 1;
                while(c)
                    if ChNum-cnt >= 1
                        A = macroMontage(ChNum-cnt).Area;
                        if strcmpi(A,ChArea)
                            cnt = cnt + 1;
                        else
                            c = 0;
                        end
                    else
                        c = 0;
                    end
                end
                ChLabel = sprintf('%s%d',ChArea,cnt);
            end
            
            
            contact_ind = [];
            for ii = 1:length(elec_name)
                if strcmpi(ChLabel(1:end-1),elec_name{ii}) && ...
                        strcmpi(ChLabel(end),num2str(elec_n(ii)))
                    contact_ind = ii;
                end
            end
            
            
            if ~isempty(contact_ind); disp(elec_name{contact_ind}); end
            if isempty(contact_ind); warning('contact missing in mloc file');
                missingContacts = [missingContacts, iContact];
                continue
            else
                elec_coord_pt_space = [x(contact_ind), y(contact_ind), z(contact_ind)];
                cfg=[];
                cfg.plotEm = 0;
                cfg.isSubdural=0; % 0 indicates that an electrode is a depth electrode
                cfg.elecCoord = elec_coord_pt_space;
                cfg.elecNames{1,1} = ChLabel;
                cfg.isLeft = strcmpi(ChLabel(1),'L');
                [avgCoords, ELEC_NAMES, isLeft]=sub2AvgBrain(anat_pt,cfg);
                
                ptID = [ptID iiP];
                PgroupAvgCoords=[PgroupAvgCoords; avgCoords];
                PgroupLabels=[PgroupLabels, ELEC_NAMES];
                PgroupIsLeft=[PgroupIsLeft; isLeft];
                contactID = [contactID, ChNum];
                PLOT_VEC = results_all(iiP).psA(iContact) < 0.05;
                NEURAL_MARKERS = [NEURAL_MARKERS, PLOT_VEC];
                avDiffCoherence = nanmean(results_all(iiP).diffCoherence(iContact,:));
                NEURAL_MARKERS_COHERENCE = [NEURAL_MARKERS_COHERENCE, avDiffCoherence];
                
                clear avgCoords ELEC_NAMES isLeft
            end
            
        end
        
    end
    
    % Display effect for targeted stimulation patients
    anatomyInfo.patients = ptID;
    anatomyInfo.NEURAL_MARKERS = NEURAL_MARKERS;
    anatomyInfo.NEURAL_MARKERS_COHERENCE = NEURAL_MARKERS_COHERENCE;
    
    anatomyInfo.PgroupAvgCoords = PgroupAvgCoords;
    anatomyInfo.PgroupLabels = PgroupLabels;
    anatomyInfo.PgroupIsLeft = PgroupIsLeft;
    anatomyInfo.behavMode = behav_str;
    
    save(fullfile(globalFsDir,'anatomyPlotsData',sprintf('mazeNeuralMarkers_%s',behav_str)),'anatomyInfo')
    
end


%% Plotting stim and probe on one brain

for ii_b = 1:length(behav_str_vec)
    behav_str = behav_str_vec{ii_b};
    anatomyInfo = load(fullfile(globalFsDir,'anatomyPlotsData',sprintf('mazeNeuralMarkers_%s',behav_str)),'anatomyInfo');
    anatomyInfo = anatomyInfo.anatomyInfo;
    NEURAL_MARKERS = anatomyInfo.NEURAL_MARKERS ;
    NEURAL_MARKERS_COHERENCE = anatomyInfo.NEURAL_MARKERS_COHERENCE ;
    PgroupAvgCoords = anatomyInfo.PgroupAvgCoords;
    PgroupLabels = anatomyInfo.PgroupLabels;
    PgroupIsLeft = anatomyInfo.PgroupIsLeft;
    ptID = anatomyInfo.patients;
    
    for ii_a = 1:2
        if (ii_a == 1)
            DIFF = NEURAL_MARKERS;
            figName = sprintf(sprintf('Maze_anatomy_Fig_ampDiff_%s',behav_str));
        else
            NEURAL_MARKERS_COHERENCE(isnan(NEURAL_MARKERS_COHERENCE)) = 0;
            DIFF = NEURAL_MARKERS_COHERENCE;
            figName = sprintf(sprintf('Maze_anatomy_Fig_CoherenceDiff_%s',behav_str));
        end
        
        f0 = figure('Name', figName,'NumberTitle','off');
        % Some WYSIWYG options:
        set(gcf,'DefaultAxesFontSize',10);
        set(gcf,'DefaultAxesFontName','arial');
        set(gcf,'PaperUnits','centimeters','PaperPosition',[0.2 0.2 19 24.7]); % this size is the maximal to fit on an A4 paper when printing to PDF
        set(gcf,'PaperOrientation','portrait');
        set(gcf,'Units','centimeters','Position', get(gcf,'paperPosition')+[1 1 0 0]);
        colormap('jet');
        
        ax1 = axes('position',[0.1,0.1,.2,.2],'units','centimeters');
        ax2 = axes('position',[0.4,0.1,.2,.2],'units','centimeters');
        ax3 = axes('position',[0.1,0.3,.2,.2],'units','centimeters');
        ax4 = axes('position',[0.4,0.3,.2,.2],'units','centimeters');
        ax5 = axes('position',[0.1,0.5,.2,.2],'units','centimeters');
        ax6 = axes('position',[0.4,0.5,.2,.2],'units','centimeters');
        ax7 = axes('position',[0.1,0.7,.2,.2],'units','centimeters');
        ax8 = axes('position',[0.4,0.7,.2,.2],'units','centimeters');
        
        N = 64;
        % map = colormap('jet');
        % map = brewermap(N,'RdYlBu');
        map = brewermap(N,'BrBG');
        map = flipud(map);
        %values = linspace(-.15,.15,N);
        values = linspace(-0.25,0.25,N);
        [colorBin,~] = discretize(DIFF,values);
        colorBin(find(DIFF < values(1))) = 1;
        colorBin(find(DIFF > values(end))) = N;
        elecColors = map(colorBin,:);
        
        
        % brain
        cfg=[];
        cfg.view='li';
        cfg.ignoreDepthElec='n';
        cfg.opaqueness=0.3;
        cfg.elecSize = 2;
        cfg.elecColors = elecColors;
        cfg.elecColorScale = [0 64];
        cfg.showLabels='n';
        cfg.title= '';
        cfg.edgeBlack = 'n';
        cfg.elecCoord = [anatomyInfo.PgroupAvgCoords,anatomyInfo.PgroupIsLeft];
        cfg.elecNames = anatomyInfo.PgroupLabels;
        cfg.axis = ax1;
        cfgOut=plotPialSurf('fsaverage',cfg);
        colorbar(ax1,'off')
        cfg.view='ri';
        cfg.axis = ax2;
        cfgOut=plotPialSurf('fsaverage',cfg);
        colorbar(ax2,'off')
        
        cfg=[];
        cfg.ignoreDepthElec='n';
        cfg.opaqueness=0.3;
        cfg.elecSize = 3;
        cfg.elecColors = elecColors;
        cfg.elecColorScale = [0 64];
        cfg.showLabels='n';
        cfg.title= '';
        cfg.edgeBlack = 'n';
        cfg.elecCoord = [anatomyInfo.PgroupAvgCoords,anatomyInfo.PgroupIsLeft];
        cfg.elecNames = anatomyInfo.PgroupLabels;
        cfg.axis = ax3;
        cfg.view='lm';
        cfgOut=plotPialSurf('fsaverage',cfg);
        cfg.view='rm';
        cfg.axis = ax4;
        cfgOut=plotPialSurf('fsaverage',cfg);
        
        
        cfg=[];
        cfg.ignoreDepthElec='n';
        cfg.opaqueness=0.3;
        cfg.elecSize = 3;
        cfg.elecColors = elecColors;
        cfg.elecColorScale = [0 64];
        cfg.showLabels='n';
        cfg.title= '';
        cfg.edgeBlack = 'n';
        cfg.elecCoord = [anatomyInfo.PgroupAvgCoords,anatomyInfo.PgroupIsLeft];
        cfg.elecNames = anatomyInfo.PgroupLabels;
        cfg.axis = ax5;
        cfg.view='lsv';
        cfgOut=plotPialSurf('fsaverage',cfg);
        cfg.view='rsv';
        cfg.axis = ax6;
        cfgOut=plotPialSurf('fsaverage',cfg);
        
        
        brainView.light=[1 0 0];
        brainView.hem='r';
        brainView.eyes=[45 0]
        cfg.view=brainView
        
        cfg=[];
        cfg.ignoreDepthElec='n';
        cfg.opaqueness=0.3;
        cfg.elecSize = 3;
        cfg.elecColors = elecColors;
        cfg.elecColorScale = [0 64];
        cfg.showLabels='n';
        cfg.title= '';
        cfg.edgeBlack = 'n';
        cfg.elecCoord = [anatomyInfo.PgroupAvgCoords,anatomyInfo.PgroupIsLeft];
        cfg.elecNames = anatomyInfo.PgroupLabels;
        cfg.axis = ax7;
        cfg.view='l';
        cfgOut=plotPialSurf('fsaverage',cfg);
        view(ax7,[-94,27])
        
        cfg.view='r';
        cfg.axis = ax8;
        cfgOut=plotPialSurf('fsaverage',cfg);
        view(ax8,[94,27])
        
        outputFigureFolder = fullfile(globalFsDir,'anatomyPlotsFigures');
        
        a = gcf;
        set(f0,'renderer','zbuffer');
        % res = 400;
        % eval(['print ', [outputFigureFolder,'\',figName], ' -f', num2str(a.Number),sprintf(' -depsc  -r%d',res), '-cmyk' ]); % adding r600 slows down this process significantly!
        res = 600;
        eval(['print ', [outputFigureFolder,'\',figName], ' -f', num2str(a.Number),sprintf(' -dtiff  -r%d',res), '-cmyk' ]); % adding r600 slows down this process significantly!
        
        figure
        figName = sprintf('rendererColorbar_%s',behav_str);
        imagesc(1*ones(1,64),1:64,values)
        colormap(map)
        cc = colorbar;
        cc.TickDirection = 'out';
        cc.Ticks = [values(1),0,values(end)];
        title(sprintf('axis limits- [%f,%f]',values(1),values(end)))
        a = gcf;
        res = 600;
        eval(['print ', [outputFigureFolder,'\',figName], ' -f', num2str(a.Number),sprintf(' -dtiff  -r%d',res), '-cmyk' ]); % adding r600 slows down this process significantly!
        
    end
    
end
%%
%% Plotting stim and probe on one brain
anatomyInfo = load(fullfile(globalFsDir,'anatomyPlotsData','stimEffectSleepOsc'),'anatomyInfo');
anatomyInfo = anatomyInfo.anatomyInfo;
ptID = anatomyInfo.patients;
anatomyInfo.patients = ptID;
allSWEvent = anatomyInfo.allSWEvent;
allSpEvent = anatomyInfo.allSpEvent;
DIFF =  allSWEvent(:,1)-allSWEvent(:,2);

figName = sprintf('stimEffects_locations_SW');
f0 = figure('Name', figName,'NumberTitle','off');
% Some WYSIWYG options:
set(gcf,'DefaultAxesFontSize',10);
set(gcf,'DefaultAxesFontName','arial');
set(gcf,'PaperUnits','centimeters','PaperPosition',[0.2 0.2 19 24.7]); % this size is the maximal to fit on an A4 paper when printing to PDF
set(gcf,'PaperOrientation','portrait');
set(gcf,'Units','centimeters','Position', get(gcf,'paperPosition')+[1 1 0 0]);
colormap('jet');

ax1 = axes('position',[0.1,0.1,.2,.2],'units','centimeters');
ax2 = axes('position',[0.4,0.1,.2,.2],'units','centimeters');
ax3 = axes('position',[0.1,0.3,.2,.2],'units','centimeters');
ax4 = axes('position',[0.4,0.3,.2,.2],'units','centimeters');
ax5 = axes('position',[0.1,0.5,.2,.2],'units','centimeters');
ax6 = axes('position',[0.4,0.5,.2,.2],'units','centimeters');
ax7 = axes('position',[0.1,0.7,.2,.2],'units','centimeters');
ax8 = axes('position',[0.4,0.7,.2,.2],'units','centimeters');

N = 64;
% map = colormap('jet');
% map = brewermap(N,'RdYlBu');
map = brewermap(N,'RdBu');
map = flipud(map);
values = linspace(-0.5,0.5,N);
[colorBin,~] = discretize(DIFF,values);
colorBin(find(DIFF < values(1))) = 1;
elecColors = map(colorBin,:);

% brain
cfg=[];
cfg.view='li';
cfg.ignoreDepthElec='n';
cfg.opaqueness=0.3;
cfg.elecSize = 2;
cfg.elecColors = elecColors;
cfg.elecColorScale = [0 64];
cfg.showLabels='n';
cfg.title= '';
cfg.edgeBlack = 'n';
cfg.elecCoord = [anatomyInfo.PgroupAvgCoords,anatomyInfo.PgroupIsLeft];
cfg.elecNames = anatomyInfo.PgroupLabels;
cfg.axis = ax1;
cfgOut=plotPialSurf('fsaverage',cfg);
colorbar(ax1,'off')
cfg.view='ri';
cfg.axis = ax2;
cfgOut=plotPialSurf('fsaverage',cfg);
colorbar(ax2,'off')

cfg=[];
cfg.ignoreDepthElec='n';
cfg.opaqueness=0.3;
cfg.elecSize = 3;
cfg.elecColors = elecColors;
cfg.elecColorScale = [0 64];
cfg.showLabels='n';
cfg.title= '';
cfg.edgeBlack = 'n';
cfg.elecCoord = [anatomyInfo.PgroupAvgCoords,anatomyInfo.PgroupIsLeft];
cfg.elecNames = anatomyInfo.PgroupLabels;
cfg.axis = ax3;
cfg.view='lm';
cfgOut=plotPialSurf('fsaverage',cfg);
cfg.view='rm';
cfg.axis = ax4;
cfgOut=plotPialSurf('fsaverage',cfg);


cfg=[];
cfg.ignoreDepthElec='n';
cfg.opaqueness=0.3;
cfg.elecSize = 3;
cfg.elecColors = elecColors;
cfg.elecColorScale = [0 64];
cfg.showLabels='n';
cfg.title= '';
cfg.edgeBlack = 'n';
cfg.elecCoord = [anatomyInfo.PgroupAvgCoords,anatomyInfo.PgroupIsLeft];
cfg.elecNames = anatomyInfo.PgroupLabels;
cfg.axis = ax5;
cfg.view='lsv';
cfgOut=plotPialSurf('fsaverage',cfg);
cfg.view='rsv';
cfg.axis = ax6;
cfgOut=plotPialSurf('fsaverage',cfg);


brainView.light=[1 0 0];
brainView.hem='r';
brainView.eyes=[45 0]
cfg.view=brainView

cfg=[];
cfg.ignoreDepthElec='n';
cfg.opaqueness=0.3;
cfg.elecSize = 3;
cfg.elecColors = elecColors;
cfg.elecColorScale = [0 64];
cfg.showLabels='n';
cfg.title= '';
cfg.edgeBlack = 'n';
cfg.elecCoord = [anatomyInfo.PgroupAvgCoords,anatomyInfo.PgroupIsLeft];
cfg.elecNames = anatomyInfo.PgroupLabels;
cfg.axis = ax7;
cfg.view='l';
cfgOut=plotPialSurf('fsaverage',cfg);
view(ax7,[-94,27])

cfg.view='r';
cfg.axis = ax8;
cfgOut=plotPialSurf('fsaverage',cfg);
view(ax8,[94,27])

outputFigureFolder = fullfile(globalFsDir,'anatomyPlotsFigures');

a = gcf;
set(f0,'renderer','zbuffer');
% res = 400;
% eval(['print ', [outputFigureFolder,'\',figName], ' -f', num2str(a.Number),sprintf(' -depsc  -r%d',res), '-cmyk' ]); % adding r600 slows down this process significantly!
res = 600;
eval(['print ', [outputFigureFolder,'\',figName], ' -f', num2str(a.Number),sprintf(' -dtiff  -r%d',res), '-cmyk' ]); % adding r600 slows down this process significantly!

figure
figName = 'rendererColorbar';
imagesc(1*ones(1,64),1:64,values)
colormap(map)
cc = colorbar;
cc.TickDirection = 'out';
cc.Ticks = [-2,0,2];
title(sprintf('axis limits- [%f,%f]',values(1),values(end)))
a = gcf;
res = 600;
eval(['print ', [outputFigureFolder,'\',figName], ' -f', num2str(a.Number),sprintf(' -dtiff  -r%d',res), '-cmyk' ]); % adding r600 slows down this process significantly!


%% Final version


anatomyInfo = load(fullfile(globalFsDir,'anatomyPlotsData','stimEffectSleepOsc'),'anatomyInfo');
anatomyInfo = anatomyInfo.anatomyInfo;
ptID = anatomyInfo.patients;
anatomyInfo.patients = ptID;
allSp_prob = anatomyInfo.allSp_prob;
DIFF =  allSp_prob(:,1)-allSp_prob(:,2);


frontalVec = logical(zeros(1,length(ptID)));
for ii = 1:length(ptID)
    area = classifyArea(anatomyInfo.PgroupLabels{ii}(5:end-1));
    if area.isFrontal
        frontalVec(ii) = true;
    end
end
[p,h,stats] = ranksum(DIFF(frontalVec),DIFF(~frontalVec));


figName = sprintf('stimEffects_locations_SW');
f0 = figure('Name', figName,'NumberTitle','off');
% Some WYSIWYG options:
set(gcf,'DefaultAxesFontSize',10);
set(gcf,'DefaultAxesFontName','arial');
set(gcf,'PaperUnits','centimeters','PaperPosition',[0.2 0.2 19 24.7]); % this size is the maximal to fit on an A4 paper when printing to PDF
set(gcf,'PaperOrientation','portrait');
set(gcf,'Units','centimeters','Position', get(gcf,'paperPosition')+[1 1 0 0]);
colormap('jet');

ax1 = axes('position',[0.1,0.1,.2,.2],'units','centimeters');
ax2 = axes('position',[0.4,0.1,.2,.2],'units','centimeters');
ax3 = axes('position',[0.1,0.3,.2,.2],'units','centimeters');
ax4 = axes('position',[0.4,0.3,.2,.2],'units','centimeters');
ax5 = axes('position',[0.1,0.5,.2,.2],'units','centimeters');
ax6 = axes('position',[0.4,0.5,.2,.2],'units','centimeters');
ax7 = axes('position',[0.1,0.7,.2,.2],'units','centimeters');
ax8 = axes('position',[0.4,0.7,.2,.2],'units','centimeters');

N = 64;
% map = colormap('jet');
% map = brewermap(N,'RdYlBu');
map = brewermap(N,'RdBu');
map = flipud(map);
values = linspace(-.1,.1,N);
[colorBin,~] = discretize(DIFF,values);
colorBin(find(DIFF < values(1))) = 1;
colorBin(find(DIFF > values(end))) = N;
elecColors = map(colorBin,:);


cfg=[];
cfg.ignoreDepthElec='n';
cfg.opaqueness=0.3;
cfg.elecSize = 3;
cfg.elecColors = elecColors;
cfg.elecColorScale = [0 64];
cfg.showLabels='n';
cfg.title= '';
cfg.edgeBlack = 'n';
cfg.elecCoord = [anatomyInfo.PgroupAvgCoords,anatomyInfo.PgroupIsLeft];
cfg.elecNames = anatomyInfo.PgroupLabels;
cfg.axis = ax7;
cfg.view='l';
cfgOut=plotPialSurf('fsaverage',cfg);
view(ax7,[-94,27])

cfg.view='r';
cfg.axis = ax8;
cfgOut=plotPialSurf('fsaverage',cfg);
view(ax8,[94,27])

outputFigureFolder = fullfile(globalFsDir,'anatomyPlotsFigures');

a = gcf;
set(f0,'renderer','zbuffer');
% res = 400;
% eval(['print ', [outputFigureFolder,'\',figName], ' -f', num2str(a.Number),sprintf(' -depsc  -r%d',res), '-cmyk' ]); % adding r600 slows down this process significantly!
res = 600;
eval(['print ', [outputFigureFolder,'\',figName], ' -f', num2str(a.Number),sprintf(' -dtiff  -r%d',res), '-cmyk' ]); % adding r600 slows down this process significantly!

figure
figName = 'rendererColorbar';
imagesc(1*ones(1,64),1:64,values)
colormap(map)
cc = colorbar;
cc.TickDirection = 'out';
cc.Ticks = [-2,0,2];
title(sprintf('axis limits- [%f,%f]',values(1),values(end)))
a = gcf;
res = 600;
eval(['print ', [outputFigureFolder,'\',figName], ' -f', num2str(a.Number),sprintf(' -dtiff  -r%d',res), '-cmyk' ]); % adding r600 slows down this process significantly!
