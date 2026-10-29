%RUN_PAPER  Reproduce every result in the MinEl_PPG paper.
%
%   1. Detection: every arm on every recording, timed end to end
%      (proposed detector, MSPTD_PB, MSPTDfast without pre-filter, and the
%      ablation arms PLAIN and PLAINREF), against the ECG or annotation
%      references. Beat detection is identical for every database: each arm
%      runs once on the whole recording, from raw PPG to final beat indices.
%      Databases differ only in references and analysis.
%   2. Per-database tables (printed), results saved to ppg_benchmark_ACC.mat.
%   3. All statistics (benchmark/run_statistics.m).
%
%   Requirements: MATLAB R2024b with Signal Processing and Statistics and
%   Machine Learning Toolboxes, the ppg-beats toolbox on the path, and the
%   databases listed in README.md. Set the file locations below.
clear; close all; clc;
addpath(fileparts(mfilename('fullpath')));
addpath(fullfile(fileparts(mfilename('fullpath')), 'benchmark'));

% Download data
download_data('.');

%% ----------------------------- CONFIG --------------------------------
cfg = struct();
cfg.DB   = {'mimic_large_test','mimic_test','mimic_af','mimic_train','bidmc','capnobase','dalia'};
cfg.DETECTORS = {'msptd_pb','msptd','user','plain','plainref'};
cfg.REFS = {'pt','minel','titen'};   % ECG references; PT is primary (column 1)
cfg.MATCH_TOL = 0.15;                % s, detection-to-reference matching tolerance

% timing: minimum of TIMING_REPS runs per recording after WARMUP_REPS warm-up runs
cfg.WARMUP_REPS = 2;
cfg.TIMING_REPS = 3;

% pulse-transit offset (windowed median; window set per database in metrics.m)
cfg.LAG_WIN_SEC          = 300;
cfg.LAG_SEARCH_INITIAL_S = 0.8;
cfg.LAG_SEARCH_LADDER_S  = [0.2 0.5 0.8];
cfg.LAG_MIN_MATCHES      = 5;

% data files
DATA = './';
cfg.files.mimic_large_test = fullfile(DATA, 'mimic-large-test/mimic_perform_large_test_a_data.mat');
cfg.files.mimic_test       = fullfile(DATA, 'mimic_test/mimic_perform_test_all_data.mat');
cfg.files.mimic_train      = fullfile(DATA, 'mimic_train/mimic_perform_train_all_data.mat');
cfg.files.mimic_af_af      = fullfile(DATA, 'mimic_perform_af/mimic_perform_af_data.mat');
cfg.files.mimic_af_nonaf   = fullfile(DATA, 'mimic_perform_af/mimic_perform_non_af_data.mat');
cfg.files.bidmc            = fullfile(DATA, 'bidmc/bidmc_data.mat');
cfg.files.capnobase        = fullfile(DATA, 'capnobase/capnobase_db.mat');
cfg.dalia.data_dir         = fullfile(DATA, 'dalia/dalia_mat');

cfg.dalia.subjects  = 1:15;
cfg.dalia.activities = {'transient','sitting','stairs','table_soccer','cycling', ...
                        'driving','lunch','walking','working'};   % IDs 0..8

M = metrics();

%% ----------------------------- DETECTION -----------------------------
for db = 1:numel(cfg.DB)
    DB = cfg.DB{db};
    [subj, Fs, is_ecg] = load_databases(DB, cfg);
    isDalia = strcmp(DB, 'dalia');
    cfg.LAG_WIN_SEC = M.lag_window(DB);
    if isDalia, refs = {'ecg_annot'}; elseif ~is_ecg, refs = {'annot'}; else, refs = cfg.REFS; end
    dets = cfg.DETECTORS; nr = numel(refs); nd = numel(dets); n = numel(subj);

    pids = {subj.pid};
    if all(cellfun(@isempty, pids)), n_pat = n; else, n_pat = numel(unique([subj.pid])); end
    fprintf('\n\n==============================================================\n');
    fprintf(' DATABASE: %s (%d patients, %d recordings, Fs = %d Hz)\n', DB, n_pat, n, Fs);
    fprintf(' Reference columns: %s\n', strjoin(refs, ' | '));
    fprintf('==============================================================\n');

    % warm-up: fills filter caches and JIT before anything is timed
    fprintf('  [warmup] %d rep(s) ... ', cfg.WARMUP_REPS);
    wp = subj(1).ppg(:); wp(~isfinite(wp)) = 0;
    nw = min(numel(wp), round(60*Fs)); wp = wp(1:nw);
    for rep = 1:cfg.WARMUP_REPS
        if is_ecg, ecg_references(subj(1).ecg(1:min(end, nw)), Fs); end
        for dt = 1:nd, run_arm(dets{dt}, wp, Fs); end
    end
    fprintf('done.\n');

    P = struct('group', cell(n,1), 'refs', cell(n,1), 'dets', cell(n,1), ...
               'cpu', cell(n,1), 'dur', cell(n,1), 'pid', cell(n,1));
    for i = 1:n
        ppg = subj(i).ppg(:); ppg(~isfinite(ppg)) = 0;
        P(i).group = subj(i).group;
        P(i).dur   = numel(ppg) / Fs;
        P(i).pid   = subj(i).pid; if isempty(P(i).pid), P(i).pid = i; end

        % references
        if is_ecg
            r = ecg_references(subj(i).ecg, Fs);
            P(i).refs = cellfun(@(f) r.(f), refs, 'UniformOutput', false);
        else
            a = subj(i).ecg(:);
            P(i).refs = {sort(unique(a(a >= 1 & a <= numel(ppg))))};
        end

        % detectors, timed: minimum of TIMING_REPS runs
        P(i).dets = cell(1, nd); P(i).cpu = zeros(1, nd);
        for dt = 1:nd
            best = Inf;
            for rep = 1:cfg.TIMING_REPS
                [pk, cpu] = run_arm(dets{dt}, ppg, Fs);
                if cpu < best, best = cpu; pk_best = pk; end
            end
            P(i).dets{dt} = pk_best(:); P(i).cpu(dt) = best;
        end
        if isDalia
            fprintf('Subj %3d: ref=[%s] det=[%s] cpu=[%s]\n', i, ...
                sprintf('%d ', cellfun(@numel, P(i).refs)), ...
                sprintf('%d ', cellfun(@numel, P(i).dets)), sprintf('%.2f ', P(i).cpu));
        end
    end

    %% ------------------------- METRICS -------------------------------
    R = struct('group', cell(n,1), 'M', cell(n,1));
    for i = 1:n
        Mi = cell(nr, nd);
        for r = 1:nr
            for dt = 1:nd
                Mi{r,dt} = M.compute(P(i).refs{r}, P(i).dets{dt}, Fs, cfg.MATCH_TOL, P(i).cpu(dt), P(i).dur, cfg);
            end
        end
        R(i).group = P(i).group; R(i).M = Mi;
    end

    gn = M.group_names(DB);
    subsets = {'ALL'}; idx_sub = {1:n};
    if ~isempty(gn) && any([R.group]==1) && any([R.group]==2)
        subsets = {'ALL', ['GROUP1 (' gn{1} ')'], ['GROUP2 (' gn{2} ')']};
        idx_sub = {1:n, find([R.group]==1), find([R.group]==2)};
    end
    for s = 1:numel(subsets)
        A = cell(nd, nr);
        for dt = 1:nd
            for r = 1:nr
                A{dt,r} = M.aggregate(arrayfun(@(k) R(k).M{r,dt}, idx_sub{s}, 'UniformOutput', false));
            end
        end
        M.print_table(sprintf('%s (%s, N=%d)', DB, subsets{s}, numel(idx_sub{s})), refs, dets, A);
    end

    % reference-agreement diagnostic
    if is_ecg && nr == 3
        c = cell2mat(arrayfun(@(p) cellfun(@numel, p.refs), P, 'UniformOutput', false));
        spread = (max(c, [], 2) - min(c, [], 2)) ./ max(c, [], 2);
        bad = find(spread > 0.10);
        fprintf('\n[%s] Recordings with >10%% reference disagreement: %d / %d\n', DB, numel(bad), n);
        for j = bad'
            fprintf(' S%3d: %s=%d %s=%d %s=%d (max diff %.0f%%)\n', j, refs{1}, c(j,1), ...
                refs{2}, c(j,2), refs{3}, c(j,3), 100*spread(j));
        end
    end

    %% ------------------------- COLLECT -------------------------------
    % reference clock: 64 Hz PPG timeline for PPG-DaLiA, the shared ECG/PPG
    % clock elsewhere (used for the non-inferiority margin)
    if ~exist('ACC', 'var')
        ACC = struct('name',{},'R',{},'P',{},'Fs',{},'refs',{},'dets',{},'tol',{},'is_ecg',{},'Fs_ref',{});
    end
    R_acc = R;
    pid = [P.pid];
    if numel(unique(pid)) < numel(pid), R_acc = collapse_by_patient(R, pid); end
    ACC(end+1) = struct('name',DB,'R',R_acc,'P',P,'Fs',Fs,'refs',{refs},'dets',{dets}, ...
                        'tol',cfg.MATCH_TOL,'is_ecg',is_ecg,'Fs_ref',Fs); %#ok<SAGROW>
    if ~isempty(gn) && any([R.group]==1) && any([R.group]==2)
        for q = 1:2
            ig = find([R.group] == q);
            ACC(end+1) = struct('name',[DB '@' gn{q}],'R',R(ig),'P',P(ig),'Fs',Fs,'refs',{refs}, ...
                'dets',{dets},'tol',cfg.MATCH_TOL,'is_ecg',is_ecg,'Fs_ref',Fs); %#ok<SAGROW>
        end
    end

    if isDalia
        ACC = dalia_activities(ACC, P, subj, Fs, refs, dets, cfg);
        durT = sum([P.dur]); cpuAll = reshape([P.cpu], nd, []);
        fprintf('\n=== DaLiA timing (%.1f h of PPG) ===\n', durT/3600);
        for dt = 1:nd
            fprintf('%-11s : %5.3f s  RTF %8.0f\n', upper(dets{dt}), sum(cpuAll(dt,:)), durT/sum(cpuAll(dt,:)));
        end
    end
end

save('ppg_benchmark_ACC.mat', 'ACC', 'cfg');
fprintf('\nResults saved to ppg_benchmark_ACC.mat\n');

%% ----------------------------- STATISTICS ----------------------------
run_statistics(ACC, cfg);

%% ----------------------------- ENVIRONMENT ---------------------------
fprintf('\nMATLAB %s\n', version);
v = ver('signal'); fprintf('Signal Processing Toolbox %s\n', v.Version);
v = ver('stats');  fprintf('Statistics and Machine Learning Toolbox %s\n', v.Version);
fprintf('Platform: %s | threads: %d | cores: %d\n', computer('arch'), maxNumCompThreads, feature('numcores'));
fprintf('All done.\n');


%% =====================================================================
%  LOCAL FUNCTIONS
%  =====================================================================
function [pk, cpu] = run_arm(arm, ppg, Fs)
% One arm on one whole recording, timed from raw PPG to final beat indices.
% Identical for every database.
switch arm
    case 'user'
        t = tic; pk = minel_ppg(ppg, Fs); cpu = toc(t);
    case 'msptd_pb'
        t = tic; pk = comparators('msptd_pb', ppg, Fs); cpu = toc(t);
    case 'msptd'
        t = tic; pk = comparators('msptd', ppg, Fs); cpu = toc(t);
    case {'plain', 'plainref'}
        opts = struct('stream', 'uniform', 'refine', strcmp(arm, 'plainref'));
        t = tic; pk = minel_ppg(ppg, Fs, opts); cpu = toc(t);
    otherwise
        error('Unknown arm: %s', arm);
end
end

function ACC = dalia_activities(ACC, P, subj, Fs, refs, dets, cfg)
% Per-activity results: beats restricted to samples of each activity.
M = metrics();
names = cfg.dalia.activities; nS = numel(P); nd = numel(dets);
act_at = @(idx, act) act(min(max(floor((idx(:) - 1) / Fs * 4) + 1, 1), numel(act)));
for a = 0:numel(names)-1
    R = struct('group', cell(nS,1), 'M', cell(nS,1)); Psub = P;
    for i = 1:nS
        act = subj(i).act;
        rr = P(i).refs{1}; rr = rr(act_at(rr, act) == a);
        dd = cellfun(@(x) x(act_at(x, act) == a), P(i).dets, 'UniformOutput', false);
        dur = sum(act == a) / 4;
        Mi = cell(1, nd);
        for dt = 1:nd
            Mi{1, dt} = M.compute(rr, dd{dt}, Fs, cfg.MATCH_TOL, P(i).cpu(dt) * dur / max(P(i).dur, eps), dur, cfg);
        end
        R(i).group = 0; R(i).M = Mi;
        Psub(i).refs = {rr}; Psub(i).dets = dd; Psub(i).dur = dur;
    end
    A = cell(nd, 1);
    for dt = 1:nd
        A{dt, 1} = M.aggregate(arrayfun(@(k) R(k).M{1, dt}, 1:nS, 'UniformOutput', false));
    end
    M.print_table(sprintf('dalia (%s, N=%d)', names{a+1}, nS), refs, dets, A);
    ACC(end+1) = struct('name', ['dalia_' names{a+1}], 'R', R, 'P', Psub, 'Fs', Fs, ...
        'refs', {refs}, 'dets', {dets}, 'tol', cfg.MATCH_TOL, 'is_ecg', false, 'Fs_ref', Fs); %#ok<AGROW>
end
end

function Rp = collapse_by_patient(R, pid)
% Held-out set: one value per patient, the mean of their recordings.
[u, ~, g] = unique(pid(:));
[nr, nd] = size(R(1).M);
Rp = struct('group', cell(numel(u),1), 'M', cell(numel(u),1));
for p = 1:numel(u)
    idx = find(g == p);
    Rp(p).group = R(idx(1)).group;
    Mi = cell(nr, nd);
    for r = 1:nr
        for d = 1:nd
            Ms = cellfun(@(k) R(k).M{r,d}, num2cell(idx), 'UniformOutput', false);
            Ms = [Ms{:}]; m = Ms(1);
            for f = fieldnames(m)'
                v = [Ms.(f{1})];
                if any(strcmp(f{1}, {'nref','ndet','cpu','dur'})), m.(f{1}) = sum(v);
                else, m.(f{1}) = mean(v, 'omitnan'); end
            end
            Mi{r,d} = m;
        end
    end
    Rp(p).M = Mi;
end
end