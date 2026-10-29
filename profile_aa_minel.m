%PROFILE_MINEL  Where does minel_ppg spend its time?
%
%   Profiles minel_ppg.m on a short-record database (MIMIC Large Testing,
%   796 recordings of about 30 s) and a long-record one (PPG-DaLiA, 15
%   recordings of about 2.4 h), and prints the 15 functions with the most
%   self time (time inside the function itself, excluding the functions it
%   calls). Local functions appear as minel_ppg>detect_beats etc.
%
%   Place in MinEl_PPG/aa_test/ (next to run_aa_test.m) and run.
clear; clc;
here = fileparts(mfilename('fullpath'));
repo = fileparts(here);
addpath(here); addpath(repo); addpath(fullfile(repo, 'benchmark'));
DATA = fullfile(repo, './');                 % same data location as run_aa_test.m

% ----------------------------- CONFIG (as in run_aa_test.m) -------------
DBS = {'mimic_large_test', 'dalia'};
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

for d = 1:numel(DBS)
    [subj, Fs] = load_databases(DBS{d}, cfg);
    X = cell(numel(subj), 1); hours = 0;
    for i = 1:numel(subj)
        x = subj(i).ppg(:); x(~isfinite(x)) = 0; X{i} = x; hours = hours + numel(x) / Fs / 3600;
    end
    minel_ppg(X{1}, Fs);                     % warm-up (filter cache, JIT)

    t = tic;                                 % plain timing, profiler off
    for i = 1:numel(X), minel_ppg(X{i}, Fs); end
    t_plain = toc(t);

    profile clear; profile on
    for i = 1:numel(X), minel_ppg(X{i}, Fs); end
    profile off
    p  = profile('info');
    FT = p.FunctionTable;
    self = zeros(numel(FT), 1);
    for k = 1:numel(FT)
        ch = FT(k).Children;
        self(k) = FT(k).TotalTime - sum([ch.TotalTime]);
    end
    T = table({FT.FunctionName}', [FT.NumCalls]', [FT.TotalTime]', self, ...
        'VariableNames', {'FunctionName', 'NumCalls', 'TotalTime', 'SelfTime'});
    T = sortrows(T, 'SelfTime', 'descend');
    T.SelfShare = round(100 * T.SelfTime / sum(T.SelfTime), 1);

    fprintf('\n==== %s: %d recordings, %.1f h, %g Hz ====\n', DBS{d}, numel(X), hours, Fs);
    fprintf('Without profiler: %.3f s total, %.3f ms per recording, %.1f ms per hour of PPG\n', ...
        t_plain, 1e3 * t_plain / numel(X), 1e3 * t_plain / hours);
    disp(T(1:min(15, height(T)), :))
end
