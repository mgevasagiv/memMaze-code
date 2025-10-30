
% ERSP was within-subject comparisons between two conditions (e.g., condition 1 vs. condition 2) for each subject. 
% This means the ERSP output is already baseline-corrected and reflects the difference between the two conditions.

% collect all the ERSPs across subjects

freqs = 1:64;
times = -100:99;

for ii_a = 1:3

    if ii_a == 1
        IDS = find(TABLE_ALL.isMTL);
    elseif ii_a == 2
        IDS = find(TABLE_ALL.isMTL);
    elseif  ii_a == 3
        IDS = find(TABLE_ALL.isMTL);
    end

    all_ersp = ersp_all(IDS,:,:);
end

% Assuming you have N subjects
N = length(IDS);
all_ersp = zeros(length(freqs), length(times), N);

for i = 1:N
    all_ersp(:, :, i) = subject_ersp{i};  % each is condition1 - condition2
    all_itc(:, :, i) = subject_itc{i};  % ITC difference: cond1 - cond2
end

mean_ersp = squeeze(mean(all_ersp, 1));


% Reshape for statcond: [freq*time, subjects]
data = reshape(all_ersp, [], N)';

% Run permutation test

nFreqs = size(mean_ersp, 1);
nTimes = size(mean_ersp, 2);

pvals = zeros(nFreqs, nTimes);

for f = 1:nFreqs
    for t = 1:nTimes
        data = squeeze(all_ersp(:, f, t));
        % statcond expects a cell array of groups
        p = statcond({data', zeros(size(data))}, 'method', 'perm', 'naccu', 2000);
        pvals(f, t) = p;
    end
end

% Reshape p-values back to [freqs � times]
pvals = reshape(pvals, length(freqs), length(times));

% apply fdr correction for all points


figure;
imagesc(times, freqs, mean_ersp);
axis xy;
colorbar;
title('Mean ERSP Difference (Cond1 - Cond2)');
hold on;
contour(times, freqs, pvals < 0.05, 1, 'k'); % significance mask




% ersp_diff: [subjects x freqs x times], e.g., 21 x 64 x 200
nSubjects = size(ersp_diff, 1);
nFreqs = size(ersp_diff, 2);
nTimes = size(ersp_diff, 3);

% Step 1: Compute t-values and uncorrected p-values
tvals = zeros(nFreqs, nTimes);
pvals = zeros(nFreqs, nTimes);

for f = 1:nFreqs
    for t = 1:nTimes
        [~, p, ~, stats] = ttest(squeeze(ersp_diff(:, f, t)));
        tvals(f, t) = stats.tstat;
        pvals(f, t) = p;
    end
end

% Step 2: Threshold p-values to create binary mask
pThresh = 0.05;
sigMask = pvals < pThresh;

% Step 3: Identify clusters
conn = bwconncomp(sigMask, 4);  % 4-connected neighborhood
clusterStats = zeros(1, conn.NumObjects);
for i = 1:conn.NumObjects
    clusterStats(i) = sum(tvals(conn.PixelIdxList{i}));
end

% Step 4: Permutation testing
nPerms = 1000;
maxClusterStats = zeros(1, nPerms);

for p = 1:nPerms
    signs = randi([0 1], nSubjects, 1) * 2 - 1;  % Random sign flip
    permData = ersp_diff .* reshape(signs, [nSubjects, 1, 1]);

    permTvals = zeros(nFreqs, nTimes);
    for f = 1:nFreqs
        for t = 1:nTimes
            [~, ~, ~, stats] = ttest(squeeze(permData(:, f, t)));
            permTvals(f, t) = stats.tstat;
        end
    end

    permMask = abs(permTvals) > tinv(1 - pThresh/2, nSubjects - 1);
    permConn = bwconncomp(permMask, 4);
    maxStat = 0;
    for i = 1:permConn.NumObjects
        stat = sum(abs(permTvals(permConn.PixelIdxList{i})));
        if stat > maxStat
            maxStat = stat;
        end
    end
    maxClusterStats(p) = maxStat;
end

% Step 5: Determine significant clusters
clusterPvals = arrayfun(@(x) mean(maxClusterStats >= abs(x)), clusterStats);
finalMask = false(nFreqs, nTimes);
for i = 1:length(clusterPvals)
    if clusterPvals(i) < 0.05
        finalMask(conn.PixelIdxList{i}) = true;
    end
end

% Step 6: Plot the final cluster mask
imagesc(finalMask);
axis xy;
xlabel('Time Points');
ylabel('Frequency Bins');
title('Significant Clusters (Cluster-Based Permutation Test)');
colorbar;

