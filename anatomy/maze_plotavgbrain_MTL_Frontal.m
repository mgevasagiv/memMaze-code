global globalFsDir;
globalFsDir='E:\MAZE\anatomy\FreesurferSubjectFolders';
addpath(genpath('C:\Users\mgeva\Documents\GitHub\iELVis\'))

patients = {'da017','da018', 'da019','da019','da022','da023','ir103'};
patients_id_anat = {'sub-017','sub-018','sub-019','sub-019_2','sub-022','sub-023','ir103'};
behav_str_vec = {'G','Gm','D','X','Xm','Test'};
summarYfilesFolder = 'E:\MAZE\data_p\populationAnalysis\macroLFP_spectralAnalysis\';

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
    summaryFile = fullfile(summarYfilesFolder, sprintf('populationSummaru_spectralAnalysis_%s.mat',behav_str));
    mm = matfile(summaryFile);
    results_all = mm.results_all;
    
    missingContacts = [];
    probe_elec_coord = cell(1,1);
    contactID = [];
    contactType = [];
    
    for iiP = 1:6
        currPt = results_all(iiP).patientName;
        disp(currPt)
        ptList{iiP} = currPt;
        results_pt = results_all(iiP);
        clear hipCh MTLCh frontCh
        AREAS = results_pt.AREAS;
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
        fr_ch_num = []; 
        MTL_ch_num = [];
        % TO DO - focus on relevant channel to plot them on average brain
        if sum(MTLCh) > 0 
            MTL_ch_num = results_pt.channels(logical(MTLCh));
        end
        if sum(frontCh) > 0 
            fr_ch_num = results_pt.channels(logical(frontCh));
        end
        if isempty(MTL_ch_num) & isempty(fr_ch_num)
            continue
        end
        channels_to_plot = [MTL_ch_num, fr_ch_num];
        channelsType = [ones(1,length(MTL_ch_num)),-ones(1,length(fr_ch_num))];
            
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
        
        for iContact = 1:length(channels_to_plot)
            ChNum = channels_to_plot(iContact);
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
                contactType = [contactType channelsType(iContact)]; % 1-MTL, -1 - frontal
                PLOT_VEC = results_all(iiP).psA(iContact) < 0.05;
%                 NEURAL_MARKERS = [NEURAL_MARKERS, PLOT_VEC];
%                 avDiffCoherence = nanmean(results_all(iiP).diffCoherence(iContact,:));
%                 NEURAL_MARKERS_COHERENCE = [NEURAL_MARKERS_COHERENCE, avDiffCoherence];
%                 
                clear avgCoords ELEC_NAMES isLeft
            end
            
        end
        
    end
    
    % Display effect for targeted stimulation patients
    anatomyInfo.patients = ptID;
    anatomyInfo.contactID = contactID;
    anatomyInfo.contactType = contactType;
    %anatomyInfo.NEURAL_MARKERS = NEURAL_MARKERS;
    %anatomyInfo.NEURAL_MARKERS_COHERENCE = NEURAL_MARKERS_COHERENCE;
    isxF = find(contactType == -1);
    anatomyInfo.FgroupAvgCoords = PgroupAvgCoords(isxF,:);
    anatomyInfo.FgroupLabels = {};
    for ii = 1:length(isxF)
        anatomyInfo.FgroupLabels{ii,1} = PgroupLabels{isxF(ii)};
    end
    anatomyInfo.FgroupIsLeft = PgroupIsLeft(isxF);
    
    isxM = find(contactType == 1);
    anatomyInfo.MgroupAvgCoords = PgroupAvgCoords(isxM,:);
    anatomyInfo.MgroupLabels = {};
    for ii = 1:length(isxM)
        anatomyInfo.MgroupLabels{ii,1} = PgroupLabels{isxM(ii)};
    end
    anatomyInfo.MgroupIsLeft = PgroupIsLeft(isxM);
    
   
    save(fullfile(globalFsDir,'anatomyPlotsData',sprintf('mazeMTLFrontLoc_%s',behav_str)),'anatomyInfo')
    
end


%% Plotting stim and probe on one brain
figName = 'MAZE_frontalMTL_contacts';
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
        


elecColors = [repmat([0.1098,0.0980,0.8431],length(anatomyInfo.FgroupLabels),1);repmat([0,0,0],length(anatomyInfo.MgroupAvgCoords),1)];

cfg=[]; 
cfg.surfType = 'pial';
cfg.view='lm';
cfg.ignoreDepthElec='n';
cfg.opaqueness=0.2;
cfg.elecSize = 2;
cfg.elecColors = elecColors;
cfg.elecColorScale = [0 64];
cfg.showLabels='n';
cfg.title= '';
cfg.elecCoord = [[anatomyInfo.FgroupAvgCoords,anatomyInfo.FgroupIsLeft];[anatomyInfo.MgroupAvgCoords,anatomyInfo.MgroupIsLeft]];  
cfg.elecNames = [anatomyInfo.FgroupLabels',anatomyInfo.MgroupLabels'];

cfg.lineWidth = 0.8;
cfg.edgeBlack = 'n';
cfg.axis = ax7;
cfgOut=plotPialSurf('fsaverage',cfg);
cfg.view='rm';
cfg.axis = ax8;
cfgOut=plotPialSurf('fsaverage',cfg);

outputFigureFolder = fullfile(globalFsDir,'anatomyPlotsFigures');

a = gcf;
set(f0,'renderer','zbuffer');
res = 400;
eval(['print ', [outputFigureFolder,'\',figName], ' -f', num2str(a.Number),sprintf(' -depsc  -r%d',res), '-cmyk' ]); % adding r600 slows down this process significantly!
res = 400;
eval(['print ', [outputFigureFolder,'\',figName], ' -f', num2str(a.Number),sprintf(' -dtiff  -r%d',res), '-cmyk' ]); % adding r600 slows down this process significantly!


