%VERIFY_FAST  Is minel_ppg_fast identical to minel_ppg, and how much faster?
%
%   For every recording of all seven databases: runs minel_ppg (current)
%   and minel_ppg_fast, checks that beats and beats_int are identical, and
%   times both (minimum of REPS runs per recording, after a warm-up).
%   Place in MinEl_PPG/aa_test/ and run.
clear; clc;
here = fileparts(mfilename('fullpath'));
repo = fileparts(here);
addpath(here); addpath(repo); addpath(fullfile(repo, 'benchmark'));
DATA = fullfile(repo, './');                 % same data location as run_aa_test.m
REPS = 3;

% ----------------------------- CONFIG (as in run_aa_test.m) -------------
DBS = {'mimic_large_test','mimic_test','mimic_af','mimic_train','bidmc','capnobase','dalia'};
cfg = struct();
cand = [dir(fullfile(pwd,  '*overlap*.mat')); ...
        dir(fullfile(here, '*overlap*.mat')); ...
        dir(fullfile(repo, '*overlap*.mat'))];
assert(~isempty(cand), 'No overlap mask (*overlap*.mat) in %s, %s or %s', pwd, here, repo);
cfg.files.overlap_mask     = fullfile(cand(1).folder, cand(1).name);
cfg.files.mimic_large_test = fullfile(DATA, 'mimic-large-test/mimic_perform_large_test_a_data.mat');
cfg.files.mimic_test       = fullfile(DATA, 'mimic_test/mimic_perform_test_all_data.mat');
cfg.files.mimic_train      = fullfile(DATA, 'mimic_train/mimic_perform_train_all_data.mat');
cfg.files.mimic_af_af      = fullfile(DATA, 'mimic_perform_af/mimic_perform_af_data.mat');
cfg.files.mimic_af_nonaf   = fullfile(DATA, 'mimic_perform_af/mimic_perform_non_af_data.mat');
cfg.files.bidmc            = fullfile(DATA, 'bidmc/bidmc_data.mat');
cfg.files.capnobase        = fullfile(DATA, 'capnobase/capnobase_db.mat');
cfg.dalia.data_dir         = fullfile(DATA, 'dalia/dalia_mat');
cfg.dalia.subjects         = 1:15;
% ------------------------------------------------------------------------

fprintf('%-17s %5s %10s %12s %12s %9s\n', 'database', 'n', 'identical', 'current (s)', 'fast (s)', 'speed-up');
tot = [0 0];
for d = 1:numel(DBS)
    [subj, Fs] = load_databases(DBS{d}, cfg);
    n = numel(subj); same = 0; T = [0 0]; bad = [];
    x0 = subj(1).ppg(:); x0(~isfinite(x0)) = 0;
    minel_ppg(x0, Fs); minel_ppg_fast(x0, Fs);           % warm-up
    for i = 1:n
        x = subj(i).ppg(:); x(~isfinite(x)) = 0;
        [b1, i1] = minel_ppg(x, Fs);
        [b2, i2] = minel_ppg_fast(x, Fs);
        if isequal(b1, b2) && isequal(i1, i2), same = same + 1; else, bad(end+1) = i; end %#ok<AGROW>
        t = [Inf Inf];
        for r = 1:REPS
            s = tic; minel_ppg(x, Fs);      t(1) = min(t(1), toc(s));
            s = tic; minel_ppg_fast(x, Fs); t(2) = min(t(2), toc(s));
        end
        T = T + t;
    end
    tot = tot + T;
    fprintf('%-17s %5d %6d/%-4d %12.3f %12.3f %8.2fx\n', DBS{d}, n, same, n, T(1), T(2), T(1) / T(2));
    if ~isempty(bad)
        fprintf('   differs in recordings: %s\n', mat2str(bad(1:min(end, 20))));
    end
end
fprintf('%-17s %5s %10s %12.3f %12.3f %8.2fx\n', 'all', '', '', tot(1), tot(2), tot(1) / tot(2));
