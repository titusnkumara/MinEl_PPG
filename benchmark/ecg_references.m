function refs = ecg_references(ecg, fs)
%ECG_REFERENCES  Reference beats from the ECG.
%   refs = ecg_references(ecg, fs) returns a struct with fields
%     pt    : Pan-Tompkins on the raw ECG (primary reference; shares no code
%             with MinEl_PPG)
%     minel : MinEl on the 8-25 Hz band-passed ECG, sparse-coded at 10 Hz
%     titen : TiTen on the same sparse stream
%   All in native sample indices. Beats closer than 0.22 s are merged
%   (272 bpm ceiling), shorter than the PPG detector's refractory period so
%   the references are never the limiting factor.
REF_REFRAC = 0.22;   % s
FS_SPARSE  = 10;     % values per second

refs.minel = []; refs.titen = []; refs.pt = [];
if isempty(ecg), return; end
ecg = ecg(:);
ecg(~isfinite(ecg)) = 0;
refrac = round(REF_REFRAC * fs);

% sparse min/max coding of the band-passed ECG
[b, a] = butter(2, [8 25]/(fs/2), 'bandpass');
ecg_filt = filtfilt(b, a, ecg);
blk = max(2, round(2*fs/FS_SPARSE));
N = floor(numel(ecg_filt)/blk)*blk;
B = reshape(ecg_filt(1:N), blk, []);
[~, imin] = min(B, [], 1); [~, imax] = max(B, [], 1);
sidx = sort([(0:size(B,2)-1)*blk + imin, (0:size(B,2)-1)*blk + imax]);
sval = ecg_filt(sidx);

% MinEl
[det, ~] = MinEl(sval, FS_SPARSE);
me = sort(sidx(det)); me = me(:);
refs.minel = me([true; diff(me) > refrac]);

% TiTen
try
    [det, ~] = TiTen(sval, FS_SPARSE);
    tt = sort(sidx(det)); tt = tt(:);
    refs.titen = tt([true; diff(tt) > refrac]);
catch
    refs.titen = [];
end

% Pan-Tompkins (pan_tompkin.m, see its header for source and license)
try
    [~, pt, ~] = pan_tompkin(ecg, fs, 0);
    pt = pt(:);
    pt = sort(pt(pt >= 1 & pt <= numel(ecg)));
catch
    pt = [];
end
refs.pt = pt([true; diff(pt) > refrac]);
end


% =========================================================================
function [qrs_i_raw, qrs_amp_raw] = MinEl(ecg, fs)
% MinEl QRS detector (optimised on MIT-BIH and QT databases).
%   ecg : band-pass filtered ECG (here: its sparse min/max stream)
%   fs  : sampling frequency (here: 10 values per second)
high_factor   = 1.3;    % threshold multiplier when RR deviation > tolerance
low_factor    = 0.7;    % threshold multiplier when RR deviation <= tolerance
rr_tol        = 0.3;    % relative RR deviation tolerance
win_len_sec   = 11;     % recent window for the percentile threshold (s)
perc          = 97;     % percentile for the baseline threshold
height_thresh = 0.104;  % fraction of percentile used as base threshold
alpha         = 0.425;  % blend factor for the hybrid signal
min_dist      = 2;      % minimum distance between beats (samples)
beta          = 0.06;   % offset factor for the first threshold
win_beat      = 12;     % beat-level moving average (samples)
win_qrs       = 3;      % QRS-level moving average (samples)

qrs_amp_raw = []; qrs_i_raw = [];
if length(ecg) < win_beat + 10, return; end
ecg = ecg(:)';
ekg_abs = abs(ecg);

% stage 1: Elgendi-style blocks of interest
ma_qrs  = conv(ekg_abs, ones(1, win_qrs)/win_qrs, 'same');
ma_beat = conv(ekg_abs, ones(1, win_beat)/win_beat, 'same');
hybrid  = alpha * ekg_abs + (1 - alpha) * ma_qrs;
THR1    = ma_beat + beta * mean(ekg_abs);
d = diff([false, hybrid > THR1, false]);
starts = find(d == 1); ends = find(d == -1) - 1;
if isempty(starts), return; end
cand_locs = zeros(size(starts)); cand_amps = zeros(size(starts));
for i = 1:length(starts)
    [cand_amps(i), idx] = max(ecg(starts(i):ends(i)));
    cand_locs(i) = starts(i) + idx - 1;
end

% global percentile threshold over the last win_len_sec seconds
window_len = min(length(ecg), round(win_len_sec * fs));
base_thresh = height_thresh * prctile(ecg(end - window_len + 1 : end), perc);
valid = cand_amps > base_thresh;
cand_locs = cand_locs(valid); cand_amps = cand_amps(valid);
if isempty(cand_locs), return; end
[cand_locs, o] = sort(cand_locs); cand_amps = cand_amps(o);

% stage 2: RR-adaptive acceptance
qrs_i_raw = zeros(1, length(cand_locs)); qrs_amp_raw = zeros(1, length(cand_locs));
nacc = 0; rr_history = zeros(1, 10); rr_cnt = 0;
for i = 1:length(cand_locs)
    t = cand_locs(i); amp = cand_amps(i);
    if nacc >= 2 && rr_cnt >= 1
        med_rr = median(rr_history(1:rr_cnt));
        rr_dev = abs(t - (qrs_i_raw(nacc) + med_rr)) / med_rr;
    else
        rr_dev = NaN;
    end
    if isnan(rr_dev) || rr_dev <= rr_tol, mult = low_factor; else, mult = high_factor; end
    if nacc >= 1 && (t - qrs_i_raw(nacc)) < min_dist, continue; end
    if amp > mult * base_thresh
        nacc = nacc + 1;
        qrs_i_raw(nacc) = t; qrs_amp_raw(nacc) = amp;
        if nacc >= 2
            new_rr = qrs_i_raw(nacc) - qrs_i_raw(nacc - 1);
            if rr_cnt < 10
                rr_cnt = rr_cnt + 1; rr_history(rr_cnt) = new_rr;
            else
                rr_history(1:9) = rr_history(2:10); rr_history(10) = new_rr;
            end
        end
    end
end
qrs_i_raw = qrs_i_raw(1:nacc); qrs_amp_raw = qrs_amp_raw(1:nacc);
end


% =========================================================================
function [qrs_i_raw, qrs_amp_raw] = TiTen(ecg, fs) %#ok<INUSD>
% TiTen QRS detector, tuned for a 10 Hz min-max downsampled ECG
% (parameters maximise F1 on MIT-BIH, QT and NST databases).
integration_window    = 2;      % 0.20 s moving average of squared ECG
min_peak_distance     = 3;      % 0.30 s between candidates
init_signal_fraction  = 0.101;
signal_smoothing      = 0.1090;
searchback_multiplier = 1.66;
threshold_fraction    = 0.197;
rr_low_bound          = 0.733;
rr_high_bound         = 1.089;
rr_start_beats        = 14;
peak_search_window    = 2;      % 0.20 s raw peak search radius
searchback_guard      = 2;      % 0.20 s guard before searchback
init_noise_fraction   = 0.5;
noise_smoothing       = 0.125;
noise_memory          = 1 - noise_smoothing;
searchback_sig_alpha  = 0.25;
searchback_sig_memory = 1 - searchback_sig_alpha;
noise_threshold_frac  = 0.5;
rr_halving_factor     = 0.5;
rr_buffer_len         = 8;
signal_memory         = 1 - signal_smoothing;

ecg = ecg(:);
integrated_signal = conv(ecg.^2, ones(1, integration_window)/integration_window, 'full');
[candidate_peaks, candidate_locs] = findpeaks_mindist(integrated_signal, min_peak_distance);
if isempty(candidate_peaks), qrs_i_raw = []; qrs_amp_raw = []; return; end
num_candidates = length(candidate_peaks);

init_len = min(20, length(integrated_signal));
signal_threshold     = max(integrated_signal(1:init_len)) * init_signal_fraction;
noise_threshold      = mean(integrated_signal(1:init_len)) * init_noise_fraction;
signal_level         = signal_threshold;
noise_level          = noise_threshold;
signal_threshold_raw = max(ecg(1:init_len)) * init_signal_fraction;
noise_threshold_raw  = mean(ecg(1:init_len)) * init_noise_fraction;
signal_level_raw     = signal_threshold_raw;
noise_level_raw      = noise_threshold_raw;

beat_counter = 0; beat_counter_raw = 0;
qrs_integrated     = zeros(1, 2*num_candidates);
qrs_raw_indices    = zeros(1, 2*num_candidates);
qrs_raw_amplitudes = zeros(1, 2*num_candidates);
selected_mean_rr = 0;

for idx = 1:num_candidates
    % RR adaptation
    if beat_counter >= rr_start_beats
        start_idx = max(1, beat_counter - rr_buffer_len);
        mean_rr_interval = mean(diff(qrs_integrated(start_idx:beat_counter)));
        current_rr = qrs_integrated(beat_counter) - qrs_integrated(beat_counter - 1);
        if current_rr <= rr_low_bound*mean_rr_interval || current_rr >= rr_high_bound*mean_rr_interval
            signal_threshold     = rr_halving_factor * signal_threshold;
            signal_threshold_raw = rr_halving_factor * signal_threshold_raw;
        else
            selected_mean_rr = mean_rr_interval;
        end
    end

    % searchback after a long pause
    if selected_mean_rr > 0 && beat_counter > 0
        expected_interval = round(searchback_multiplier * selected_mean_rr);
        if (candidate_locs(idx) - qrs_integrated(beat_counter)) >= expected_interval
            search_start = qrs_integrated(beat_counter) + searchback_guard;
            search_end   = candidate_locs(idx) - searchback_guard;
            if search_start < search_end && search_start >= 1 && search_end <= length(integrated_signal)
                [peak_val, peak_idx] = max(integrated_signal(search_start:search_end));
                peak_idx = search_start + peak_idx - 1;
                if peak_val > noise_threshold
                    beat_counter = beat_counter + 1;
                    qrs_integrated(beat_counter) = peak_idx;
                    left_raw = max(1, peak_idx - peak_search_window);
                    right_raw = min(length(ecg), peak_idx);
                    [amp_raw, idx_raw] = max(ecg(left_raw:right_raw));
                    idx_raw = left_raw + idx_raw - 1;
                    if amp_raw > noise_threshold_raw
                        beat_counter_raw = beat_counter_raw + 1;
                        qrs_raw_indices(beat_counter_raw) = idx_raw;
                        qrs_raw_amplitudes(beat_counter_raw) = amp_raw;
                        signal_level_raw = searchback_sig_alpha*amp_raw + searchback_sig_memory*signal_level_raw;
                    end
                    signal_level = searchback_sig_alpha*peak_val + searchback_sig_memory*signal_level;
                end
            end
        end
    end

    % normal detection
    left_raw  = max(1, candidate_locs(idx) - peak_search_window);
    right_raw = min(length(ecg), candidate_locs(idx));
    if candidate_peaks(idx) >= signal_threshold
        [amp_raw, idx_raw] = max(ecg(left_raw:right_raw));
        idx_raw = left_raw + idx_raw - 1;
        skip_beat = false;
        if beat_counter >= 3                  % reject unusually wide complexes
            backward_width = 1;
            for k = candidate_locs(idx)-1 : -1 : max(1, candidate_locs(idx)-10)
                if integrated_signal(k) < 0.5*signal_threshold, break; end
                backward_width = backward_width + 1;
            end
            forward_width = 1;
            for k = candidate_locs(idx)+1 : min(length(integrated_signal), candidate_locs(idx)+10)
                if integrated_signal(k) < 0.5*signal_threshold, break; end
                forward_width = forward_width + 1;
            end
            if backward_width > 3 && forward_width > 3
                skip_beat = true;
                noise_level_raw = noise_smoothing*amp_raw + noise_memory*noise_level_raw;
                noise_level     = noise_smoothing*candidate_peaks(idx) + noise_memory*noise_level;
            end
        end
        if ~skip_beat
            beat_counter = beat_counter + 1;
            qrs_integrated(beat_counter) = candidate_locs(idx);
            if amp_raw >= signal_threshold_raw
                beat_counter_raw = beat_counter_raw + 1;
                qrs_raw_indices(beat_counter_raw) = idx_raw;
                qrs_raw_amplitudes(beat_counter_raw) = amp_raw;
                signal_level_raw = signal_smoothing*amp_raw + signal_memory*signal_level_raw;
            end
            signal_level = signal_smoothing*candidate_peaks(idx) + signal_memory*signal_level;
        end
    else                                      % noise (between or below thresholds)
        [amp_raw, ~] = max(ecg(left_raw:right_raw));
        noise_level_raw = noise_smoothing*amp_raw + noise_memory*noise_level_raw;
        noise_level     = noise_smoothing*candidate_peaks(idx) + noise_memory*noise_level;
    end

    % adaptive thresholds
    if noise_level ~= 0 || signal_level ~= 0
        signal_threshold = noise_level + threshold_fraction*(signal_level - noise_level);
        noise_threshold  = noise_threshold_frac * signal_threshold;
    end
    if noise_level_raw ~= 0 || signal_level_raw ~= 0
        signal_threshold_raw = noise_level_raw + threshold_fraction*(signal_level_raw - noise_level_raw);
        noise_threshold_raw  = noise_threshold_frac * signal_threshold_raw;
    end
end
qrs_i_raw   = qrs_raw_indices(1:beat_counter_raw);
qrs_amp_raw = qrs_raw_amplitudes(1:beat_counter_raw);
end

% =========================================================================
function [pks, locs] = findpeaks_mindist(x, minDist)
% Local maxima with a minimum separation (greedy, highest first);
% equivalent to findpeaks(x, 'MinPeakDistance', minDist).
x = x(:);
isMax = (x(2:end-1) > x(1:end-2)) & (x(2:end-1) > x(3:end));
allLocs = find(isMax) + 1;
allPks  = x(allLocs);
if isempty(allLocs), pks = []; locs = []; return; end
[~, o] = sort(allPks, 'descend');
sLocs = allLocs(o); sPks = allPks(o);
keep = true(size(sLocs));
for i = 1:length(sLocs)
    if keep(i)
        tooClose = abs(sLocs - sLocs(i)) < minDist;
        tooClose(1:i) = false;
        keep = keep & ~tooClose;
    end
end
[locs, o] = sort(sLocs(keep));
pks = sPks(keep); pks = pks(o);
end
