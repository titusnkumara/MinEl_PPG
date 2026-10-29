%COMPARE_WITH_MATLAB  Run the C detectors on real recordings and compare with MATLAB.
%
%   Build first (in this folder, from a terminal):  make
%   Then set DATA below to your data folder and run this script in MATLAB.
%
%   For each recording, the PPG is written to a binary file, the C program
%   run_detector is called, and its beats are compared with minel_ppg.m
%   (and, if ppg-beats is on the path, with MSPTD_PB from comparators.m).
%   Exact agreement is expected; any difference is listed.
%   The databases cover 64 Hz (non-integer block length), 125 Hz and
%   300 Hz (decimation to 100 Hz).
here = fileparts(mfilename('fullpath'));
addpath(fileparts(here)); addpath(fullfile(fileparts(here), 'benchmark'));
exe = fullfile(here, 'run_detector');
if ispc, exe = [exe '.exe']; end
assert(exist(exe, 'file') == 2, 'Build the C code first: run make in %s', here);

% ------------------------------ settings ------------------------------
DATA  = '../';                                    % folder holding the databases
DBS   = {'mimic_test', 'capnobase', 'dalia'};    % any of: mimic_large_test, mimic_test,
                                                 % mimic_train, mimic_af, bidmc, capnobase, dalia
N_MAX = Inf;                                     % recordings per database (Inf = all, as in the paper)
DO_PB = exist('detect_ppg_beats', 'file') == 2;  % also check MSPTD_PB if ppg-beats is available
DO_TIME = true;                                  % also time the three C detectors on each recording
% Timing uses the paper's workload: every recording processed whole, on
% every database (as run_paper.m).

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

fin = [tempname '.bin']; fout = [tempname '.bin'];
summary = {}; timing = {};
for d = 1:numel(DBS)
    DB = DBS{d};
    [subj, Fs] = load_databases(DB, cfg);
    if isempty(subj), fprintf('%s: no recordings found, skipped.\n', DB); continue; end
    n = min(N_MAX, numel(subj));
    fprintf('\n%s (%g Hz): checking %d recordings ...\n', DB, Fs, n);
    n_same = [0 0]; T = zeros(1, 3); dur = 0;
    for i = 1:n
        x = subj(i).ppg(:); x(~isfinite(x)) = 0;
        write_signal(fin, x, Fs);
        % MinEl_PPG
        run_c(exe, 'minel', fin, fout);
        bc = read_beats(fout); bm = minel_ppg(x, Fs);
        same = isequal(bc(:), bm(:)); n_same(1) = n_same(1) + same;
        if ~same, fprintf('  MinEl_PPG rec %3d: %d beats differ (C %d, MATLAB %d)\n', i, numel(setxor(bc, bm)), numel(bc), numel(bm)); end
        % MSPTD_PB
        if DO_PB
            run_c(exe, 'msptd_pb', fin, fout);
            bc = read_beats(fout); bm = comparators('msptd_pb', x, Fs);
            same = isequal(bc(:), bm(:)); n_same(2) = n_same(2) + same;
            if ~same, fprintf('  MSPTD_PB  rec %3d: %d beats differ (C %d, MATLAB %d)\n', i, numel(setxor(bc, bm)), numel(bc), numel(bm)); end
        end
        % timing (C, minimum of repeated runs), whole recording as in the paper
        if DO_TIME
            write_signal(fin, x, Fs);
            [st, out] = system(sprintf('"%s" time "%s" "%s"', exe, fin, fout));
            if st ~= 0
                error(['run_detector time failed (status %d): %s\n' ...
                       'If this mentions usage, rebuild: run make in %s'], st, strtrim(out), here);
            end
            fid = fopen(fout); t = fread(fid, 3, 'double')'; fclose(fid);   % s: MinEl_PPG, MSPTD_PB, MSPTD
            T = T + t;
            dur = dur + numel(x) / Fs;
        end
    end
    if DO_TIME
        timing{end+1} = sprintf('%-17s %5g Hz  %6.2f h | MinEl_PPG %8.0fx | MSPTD_PB %7.0fx (%5.1fx slower) | MSPTD %8.0fx (%4.2fx)', ...
            DB, Fs, dur/3600, dur/T(1), dur/T(2), T(2)/T(1), dur/T(3), T(3)/T(1)); %#ok<AGROW>
    end
    line = sprintf('%-17s %5g Hz  %3d recordings  MinEl_PPG identical %3d', DB, Fs, n, n_same(1));
    if DO_PB, line = sprintf('%s  MSPTD_PB identical %3d', line, n_same(2)); end
    summary{end+1} = line; %#ok<AGROW>
end
fprintf('\n=== C vs MATLAB ===\n%s\n', strjoin(summary, newline));
if DO_TIME
    fprintf(['\n=== C timing, real time factor (one core; minimum of repeated runs per recording) ===\n' ...
             '    MSPTD = MSPTDfast without pre-filter; (x) = its time relative to MinEl_PPG\n%s\n'], strjoin(timing, newline));
end
if exist(fin, 'file'), delete(fin); end
if exist(fout, 'file'), delete(fout); end


function write_signal(f, x, fs)
% int32 n, float64 fs, float64 x[n]
fid = fopen(f, 'w'); fwrite(fid, numel(x), 'int32'); fwrite(fid, fs, 'double'); fwrite(fid, x, 'double'); fclose(fid);
end

function b = read_beats(f)
% int32 m, int32 beats[m]
fid = fopen(f); m = fread(fid, 1, 'int32'); b = fread(fid, m, 'int32'); fclose(fid);
end

function run_c(exe, det, fin, fout)
status = system(sprintf('"%s" %s "%s" "%s"', exe, det, fin, fout));
if status ~= 0, error('run_detector failed (status %d)', status); end
end
