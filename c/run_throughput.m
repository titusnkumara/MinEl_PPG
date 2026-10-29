%RUN_THROUGHPUT  Whole-machine throughput of the C detectors on the paper's data.
%
%   Build first (in this folder, from a terminal):  make
%   Then set DATA below and run this script in MATLAB.
%
%   Every recording of the selected databases is written to a binary pack
%   (NaN set to 0, as in run_paper.m) and processed whole by the C program
%   throughput, which hands recordings to all cores from a shared queue.
%   For each thread count the whole set is processed REPS times and the
%   fastest wall-clock time is kept. Detection is identical to the paper's;
%   only the number of recordings processed at once changes.
%
%   Output: the aggregate real-time factor (seconds of PPG per second of
%   wall-clock time) at each thread count, and what it means in people:
%   one person produces one second of PPG per second, so the real-time
%   factor is the number of people whose PPG one machine keeps up with.
%   Timing covers detection with signals in memory: no file reading, no
%   network, no streaming latency.
here = fileparts(mfilename('fullpath'));
addpath(fileparts(here)); addpath(fullfile(fileparts(here), 'benchmark'));
exe = fullfile(here, 'throughput');
assert(exist(exe, 'file') == 2, 'Build the C code first: run make in %s', here);

% ------------------------------ settings ------------------------------
DATA    = '../';                                  % folder holding the databases
DBS     = {'mimic_large_test','mimic_test','mimic_af','mimic_train','bidmc','capnobase','dalia'};
DETS    = {'minel', 'msptd_pb', 'msptd', 'msptd_bp'};   % msptd = MSPTDfast without pre-filter;
                                                 % msptd_bp = MSPTDfast after MinEl_PPG's band-pass
THREADS = [1 2 4 8 10];                          % M1 Max: 8 performance + 2 efficiency cores
REPS    = 3;                                     % best of REPS passes over all recordings
FS_MIN  = 30;                                    % MinEl_PPG decimation floor (Hz, >= 30; as opts.fs_min)
PER_DB  = true;                                  % also report each database at max threads
POP     = 340e6;                                 % population to monitor (US, approx.; update as needed)
POPNAME = 'US population';

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

%% write one pack per database
packs = {}; names = {}; nrec = 0;
for d = 1:numel(DBS)
    [subj, Fs] = load_databases(DBS{d}, cfg);
    if isempty(subj), fprintf('%s: no recordings found, skipped.\n', DBS{d}); continue; end
    f = fullfile(tempdir, sprintf('ppgpack_%s.bin', DBS{d}));
    fid = fopen(f, 'w');
    fwrite(fid, numel(subj), 'int32');
    for i = 1:numel(subj)
        x = subj(i).ppg(:); x(~isfinite(x)) = 0;
        fwrite(fid, Fs, 'double'); fwrite(fid, numel(x), 'int32'); fwrite(fid, x, 'double');
    end
    fclose(fid);
    packs{end+1} = f; names{end+1} = DBS{d}; nrec = nrec + numel(subj); %#ok<AGROW>
    fprintf('%-17s %4d recordings packed\n', DBS{d}, numel(subj));
end
allpacks = strjoin(cellfun(@(p) ['"' p '"'], packs, 'UniformOutput', false), ' ');
fout = [tempname '.bin'];

%% all databases together, every thread count
res = struct();
for k = 1:numel(DETS)
    fprintf('\n');
    tlist = strjoin(arrayfun(@num2str, THREADS, 'UniformOutput', false), ',');
    cmd = sprintf('REPS=%d "%s" -f %g %s %s "%s" %s', REPS, exe, FS_MIN, DETS{k}, tlist, fout, allpacks);
    st = system(cmd);
    if st ~= 0, error('throughput failed (status %d). Rebuild with make in %s', st, here); end
    [S, T, W] = read_result(fout);
    res.(DETS{k}) = struct('S', S, 'T', T, 'W', W);
end

%% each database alone, at the largest thread count
if PER_DB
    fprintf('\n=== Each database alone, %d threads (real-time factor) ===\n', max(THREADS));
    rtf_db = zeros(numel(packs), numel(DETS));
    fprintf('%-17s %8s', 'database', 'hours'); fprintf(' %14s', DETS{:}); fprintf('\n');
    for d = 1:numel(packs)
        line = ''; Sd = 0;
        for k = 1:numel(DETS)
            [~, ~] = system(sprintf('REPS=%d "%s" -f %g %s %d "%s" "%s"', REPS, exe, FS_MIN, DETS{k}, max(THREADS), fout, packs{d}));
            [Sd, ~, Wd] = read_result(fout);
            rtf_db(d, k) = Sd / Wd;
            line = [line sprintf(' %13.0fx', Sd / Wd)]; %#ok<AGROW>
        end
        fprintf('%-17s %8.1f%s\n', names{d}, Sd / 3600, line);
    end
    fprintf('\n=== Machines like this one needed to monitor the %s (%s people) in real time, by database ===\n', POPNAME, grp(POP));
    fprintf('%-17s', 'database'); fprintf(' %14s', DETS{:}); fprintf('\n');
    for d = 1:numel(packs)
        c = arrayfun(@(v) grp(ceil(POP / v)), rtf_db(d, :), 'UniformOutput', false);
        fprintf('%-17s', names{d}); fprintf(' %14s', c{:}); fprintf('\n');
    end
end

%% familiar numbers
fprintf('\n=== What one machine can do (all %d databases, %d threads) ===\n', numel(packs), max(THREADS));
day = 86400; year = 365.25 * day;
for k = 1:numel(DETS)
    r = res.(DETS{k}); [~, j] = max(r.T); rtf = r.S / r.W(j); j1 = find(r.T == 1, 1);
    fprintf('\n%s\n', label(DETS{k}));
    fprintf('  Paper data (%.0f h of PPG, %d recordings) processed in %.2f s\n', r.S / 3600, nrec, r.W(j));
    fprintf('  Real-time factor: %s (one core: %s; speed-up %.1fx on %d threads)\n', ...
        grp(rtf), grp(r.S / r.W(j1)), r.W(j1) / r.W(j), r.T(j));
    fprintf('  People monitored continuously in real time: %s\n', grp(floor(rtf)));
    fprintf('  Machines like this one for the %s (%s people): %s\n', POPNAME, grp(POP), grp(ceil(POP / rtf)));
    fprintf('  One day of PPG from 1 million people: %s\n', dur_str(1e6 * day / rtf));
    fprintf('  One year of PPG from one person: %s\n', dur_str(year / rtf));
    fprintf('  Person-years of PPG processed per hour: %s\n', grp(round(rtf * 3600 / year)));
end
if isfield(res, 'minel') && isfield(res, 'msptd_pb')
    [~, j] = max(res.minel.T);
    fprintf('\nMinEl_PPG vs MSPTD_PB on the whole machine: %.1fx faster\n', res.msptd_pb.W(j) / res.minel.W(j));
    fprintf('\n=== Machines like this one for the %s (all databases pooled) ===\n', POPNAME);
    for k = 1:numel(DETS)
        r = res.(DETS{k}); [~, jj] = max(r.T);
        fprintf('  %-45s %10s machines\n', label(DETS{k}), grp(ceil(POP / (r.S / r.W(jj)))));
    end
end
for d = 1:numel(packs), delete(packs{d}); end
if exist(fout, 'file'), delete(fout); end


function [S, T, W] = read_result(f)
% int32 k, float64 signal seconds, k x (float64 threads, float64 wall seconds)
fid = fopen(f); k = fread(fid, 1, 'int32'); S = fread(fid, 1, 'double');
v = fread(fid, [2 k], 'double'); fclose(fid); T = v(1, :); W = v(2, :);
end

function s = label(det)
switch det
    case 'minel',    s = 'MinEl_PPG';
    case 'msptd_pb', s = 'MSPTD_PB (MSPTDfast as run by ppg-beats)';
    case 'msptd_bp', s = 'MSPTDfast with MinEl_PPG''s band-pass';
    otherwise,       s = 'MSPTDfast without pre-filter';
end
end

function s = grp(v)
% thousands separators
s = sprintf('%.0f', v); n = numel(s);
for p = n-3:-3:1, s = [s(1:p) ',' s(p+1:end)]; end
end

function s = dur_str(sec)
if sec < 120, s = sprintf('%.1f s', sec);
elseif sec < 7200, s = sprintf('%.1f min', sec / 60);
elseif sec < 2 * 86400, s = sprintf('%.1f h', sec / 3600);
else, s = sprintf('%.1f days', sec / 86400);
end
end
