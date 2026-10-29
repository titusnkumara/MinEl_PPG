%RUN_AA_TEST  Does anti-aliasing the decimation step change MinEl_PPG's results?
%
%   Runs three arms of MinEl_PPG on all seven databases, identical except
%   for how the input is reduced to the 30 Hz floor:
%     stride  keep every D-th sample (current minel_ppg.m, no anti-aliasing)
%     box1    block mean (boxcar, CIC order 1)
%     box2    two cascaded boxcars (triangle, CIC order 2)
%   Every arm processes each recording whole and is scored against the
%   primary reference (Pan-Tompkins, or the annotations on CapnoBase and
%   PPG-DaLiA) with the same matching and alignment as run_paper.m.
%
%   Output: the standard per-database tables (metrics.print_table), then
%   paired comparisons of box1 and box2 against stride (median paired
%   difference and Wilcoxon signed-rank p) for F1, sensitivity and RR RMSE;
%   on the held-out set at patient level, as in the paper.
%
%   Place this folder inside the repository (MinEl_PPG/aa_test/) and run it
%   from anywhere. Data are read from the repository's parent folder, as in
%   run_paper.m (DATA = '../').
clear; clc;
here = fileparts(mfilename('fullpath'));
repo = fileparts(here);
addpath(here); addpath(repo); addpath(fullfile(repo, 'benchmark'));
DATA = fullfile(repo, './');

%% ----------------------------- CONFIG --------------------------------
DBS    = {'mimic_large_test','mimic_test','mimic_af','mimic_train','bidmc','capnobase','dalia'};
ARMS   = {'stride', 'box1', 'box2'};
AA     = [0 1 2];               % opts.aa for each arm
FS_MIN = 30;                    % decimation floor (Hz), as in the paper
REPS   = 3;                     % timing: minimum of REPS runs per recording

cfg = struct();
cfg.MATCH_TOL = 0.15;
cfg.LAG_WIN_SEC          = 300;
cfg.LAG_SEARCH_INITIAL_S = 0.8;
cfg.LAG_SEARCH_LADDER_S  = [0.2 0.5 0.8];
cfg.LAG_MIN_MATCHES      = 5;
cfg.files.overlap_mask = fullfile(repo, 'overlap_mask.mat');   % use the file name from run_paper.m
cfg.files.mimic_large_test = fullfile(DATA, 'mimic-large-test/mimic_perform_large_test_a_data.mat');
cfg.files.mimic_test       = fullfile(DATA, 'mimic_test/mimic_perform_test_all_data.mat');
cfg.files.mimic_train      = fullfile(DATA, 'mimic_train/mimic_perform_train_all_data.mat');
cfg.files.mimic_af_af      = fullfile(DATA, 'mimic_perform_af/mimic_perform_af_data.mat');
cfg.files.mimic_af_nonaf   = fullfile(DATA, 'mimic_perform_af/mimic_perform_non_af_data.mat');
cfg.files.bidmc            = fullfile(DATA, 'bidmc/bidmc_data.mat');
cfg.files.capnobase        = fullfile(DATA, 'capnobase/capnobase_db.mat');
cfg.dalia.data_dir         = fullfile(DATA, 'dalia/dalia_mat');
cfg.dalia.subjects         = 1:15;
% ----------------------------------------------------------------------

M = metrics();
na = numel(ARMS);
summary = {};
for d = 1:numel(DBS)
    DB = DBS{d};
    [subj, Fs, is_ecg] = load_databases(DB, cfg);
    if isempty(subj), fprintf('%s: no recordings, skipped\n', DB); continue; end
    cfg.LAG_WIN_SEC = M.lag_window(DB);
    n = numel(subj);
    D = max(1, floor(Fs / FS_MIN));
    fprintf('\n==== %s: %d recordings, %g Hz -> %.2f Hz (D = %d) ====\n', DB, n, Fs, Fs / D, D);

    % warm-up (filter-design cache, JIT)
    wp = subj(1).ppg(:); wp(~isfinite(wp)) = 0; wp = wp(1:min(end, round(60 * Fs)));
    for a = 1:na, minel_ppg_aa(wp, Fs, struct('fs_min', FS_MIN, 'aa', AA(a))); end

    R = cell(n, na); pid = zeros(n, 1); grp = zeros(n, 1);
    for i = 1:n
        ppg = subj(i).ppg(:); ppg(~isfinite(ppg)) = 0;
        dur = numel(ppg) / Fs;
        if is_ecg
            r = ecg_references(subj(i).ecg, Fs); ref = r.pt;
        else
            q = subj(i).ecg(:); ref = sort(unique(q(q >= 1 & q <= numel(ppg))));
        end
        if isfield(subj, 'pid') && ~isempty(subj(i).pid), pid(i) = subj(i).pid; else, pid(i) = i; end
        if isfield(subj, 'group') && ~isempty(subj(i).group), grp(i) = subj(i).group; end
        for a = 1:na
            o = struct('fs_min', FS_MIN, 'aa', AA(a)); best = Inf;
            for rep = 1:REPS
                t = tic; pk = minel_ppg_aa(ppg, Fs, o); c = toc(t);
                if c < best, best = c; pk_best = pk; end
            end
            R{i, a} = M.compute(ref, pk_best(:), Fs, cfg.MATCH_TOL, best, dur, cfg);
        end
    end

    % standard table (all recordings)
    A = cell(na, 1);
    for a = 1:na, A{a, 1} = M.aggregate(R(:, a)'); end
    M.print_table(sprintf('%s (ALL, N=%d)', DB, n), {'ref'}, upper(ARMS), A);

    % paired comparisons against stride
    F1 = getm(R, {'f1'}); SE = getm(R, {'sens', 'se', 'sensitivity'}); RR = getm(R, {'rr_rmse', 'rrrmse', 'rr_rmse_s'});
    subsets = {'all', 1:n};
    gn = M.group_names(DB);
    if ~isempty(gn) && any(grp == 1) && any(grp == 2)
        subsets = [subsets; {gn{1}, find(grp == 1)'}; {gn{2}, find(grp == 2)'}];
    end
    for s = 1:size(subsets, 1)
        idx = subsets{s, 2};
        if numel(unique(pid(idx))) < numel(idx)          % held-out set: patient level
            [u, ~, g] = unique(pid(idx));
            col = @(X) cell2mat(arrayfun(@(k) mean(X(idx(g == k), :), 1, 'omitnan'), (1:numel(u))', 'UniformOutput', false));
            lvl = sprintf('%d patients', numel(u));
        else
            col = @(X) X(idx, :); lvl = sprintf('%d recordings', numel(idx));
        end
        f = col(F1); se = col(SE); rr = col(RR);
        fprintf('\n  %s / %s (%s)\n', DB, subsets{s, 1}, lvl);
        fprintf('  %-10s %8s %8s %8s | %-28s %-28s\n', 'metric', 'stride', 'box1', 'box2', 'box1 - stride (p)', 'box2 - stride (p)');
        show('F1', f, 3); show('Sens', se, 3); show('RR RMSE ms', 1e3 * rr, 1);
        summary{end+1} = sprintf('%-17s %-8s F1 %.3f / %.3f / %.3f   (box1 %+.4f, p=%.2g; box2 %+.4f, p=%.2g)', ...
            DB, subsets{s, 1}, median(f(:,1), 'omitnan'), median(f(:,2), 'omitnan'), median(f(:,3), 'omitnan'), ...
            pdiff(f, 2), sr(f, 2), pdiff(f, 3), sr(f, 3)); %#ok<SAGROW>
    end
    cpu = cellfun(@(m) m.cpu, R); dur_all = sum(cellfun(@(m) m.dur, R(:, 1)));
    fprintf('  real-time factor: stride %.0f, box1 %.0f, box2 %.0f\n', dur_all ./ sum(cpu, 1));
end

fprintf('\n===== SUMMARY: median F1 stride / box1 / box2 (median paired difference, Wilcoxon p) =====\n');
fprintf('%s\n', summary{:});


% =========================================================================
function X = getm(R, names)
% n x arms matrix of one per-recording metric, found by field name.
f = fieldnames(R{1, 1}); k = find(ismember(lower(f), names), 1);
if isempty(k)
    error('run_aa_test:field', 'No field named %s in metrics.compute output. Available: %s', ...
        strjoin(names, '/'), strjoin(f', ', '));
end
X = cellfun(@(m) double(m.(f{k})), R);
end

function d = pdiff(X, a)
d = median(X(:, a) - X(:, 1), 'omitnan');
end

function p = sr(X, a)
ok = all(isfinite(X(:, [1 a])), 2);
if nnz(ok) < 2 || all(X(ok, a) == X(ok, 1)), p = 1; return; end
p = signrank(X(ok, a), X(ok, 1));
end

function show(name, X, dec)
fmt = sprintf('%%8.%df', dec);
fprintf(['  %-10s ' fmt ' ' fmt ' ' fmt ' | %+9.4f (p=%-8.2g)       %+9.4f (p=%-8.2g)\n'], name, ...
    median(X(:,1), 'omitnan'), median(X(:,2), 'omitnan'), median(X(:,3), 'omitnan'), ...
    pdiff(X, 2), sr(X, 2), pdiff(X, 3), sr(X, 3));
end
