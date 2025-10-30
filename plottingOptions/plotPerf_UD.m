% Maze healthy cohort performance 

filePath = 'C:\Users\mgeva\Box\MAZE_project\DATA\BEHAVIOR\UD\hc_data_summary.csv';
T = readtable(filePath);
rows = 1:length(T.subj);
for jj = 1:length(rows)
    id(jj) = str2num(T.subj{rows(jj)}(3:end));
end
subjects = unique(id);

clear DATA
for ii = 1:3
    for jj = 1:length(subjects)
        rows = find(logical(T.repetition == (ii-1))' & logical(id == subjects(jj)));
        DATA(ii,jj) = 100*sum(T.quiz_dist(rows)==0)/length(rows);
    end
end

newA4figure('mazePerf')
cmap = brewermap(3,'Dark2');
axes('position',[0.1 0.1 0.3 0.4])
hold on
plot(1:3, DATA,'.-','color',[0.4 0.4 0.4])

for ii = 1:3
    x = DATA(ii,:)';
    my_boxplot(x,ii,cmap(ii,:))
    % violinplot(ii, x,'color',cmap(ii,:),'range',[0 max(max(x),100)])
end
plot(get(gca,'xlim'),[0 0],'k')
set(gca,'xtick',[1:3],'XTickLabel','','ylim',[0 100],'YTick',[0:50:100],'FontSize',16)
hold on

plotsFolder = 'E:\Dropbox\RanganathLab\grantWriting\image';
PrintActiveFigs(plotsFolder)


