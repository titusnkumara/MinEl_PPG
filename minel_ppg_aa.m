function [beats, beats_int] = minel_ppg_aa(ppg, fs, opts)
%MINEL_PPG_AA  minel_ppg with optional block-average anti-aliasing (test version).
%
%   Identical to minel_ppg.m except step 1, which can average each block of
%   D samples (boxcar, CIC order P) instead of keeping every D-th sample:
%     opts.aa = 0   keep every D-th sample (identical to minel_ppg.m)
%     opts.aa = 1   block mean (boxcar): nulls at multiples of fs/D
%     opts.aa = 2   two cascaded boxcars (triangle): twice the attenuation in dB
%   The filter's delay, aa*(D-1)/2 native samples, is removed when beats are
%   mapped back to native indices.
%
%   --- original help follows ---
%MINEL_PPG  PPG beat detection from sparse extremum coding.
%
%   beats = minel_ppg(ppg, fs)
%
%   ppg   : raw PPG signal (vector)
%   fs    : sampling frequency in Hz
%   beats : beat locations as native sample indices (column vector),
%           at the point of steepest upslope of each pulse
%
%   Pipeline (only steps 1-3 run at the native sampling rate):
%     1. The signal is reduced by keeping every D-th sample, where D is the
%        largest integer that leaves at least fs_min (default 30 Hz); the
%        reduced rate fs/D need not be an integer (125 Hz -> 31.25 Hz).
%     2. The signal is differenced and band-passed (0.5-5 Hz, first-order
%        zero-phase Butterworth).
%     3. Sparse extremum coding: the derivative is split into 0.2 s blocks and
%        only each block's minimum and maximum are kept, with their sample
%        indices (ten values per second at any sampling rate).
%     4. Rhythm-aware detection on the sparse stream (detect_beats).
%     5. Adaptive refractory blanking (blank_refractory).
%     6. Slope and amplitude quality gate (quality_gate).
%     7. Three-point parabolic refinement, mapped to the nearest native sample.
%
%   Steps 4-6 work on positions in the sparse stream; sample indices are
%   looked up from sidx only where needed.
%
%   [beats, beats_int] = minel_ppg(...) also returns the beats before step 7
%   (integer indices of the sparse extrema, in native samples).
%
%   minel_ppg(ppg, fs, opts) options:
%     opts.fs_min  lowest rate kept by step 1, in Hz (default 30, minimum 30).
%                  A higher value keeps more samples, e.g. 100.
%     opts.stream = 'uniform' runs the same detector on a uniformly resampled
%                  10 Hz stream instead of the sparse one (ablation arm PLAIN);
%                  opts.refine = true adds the +/-50 ms native-rate refinement
%                  (PLAINREF). These are for the ablation study only.
%
%   Reference: MinEl_PPG: Sparse Min-Max Coding for Fast, Accurate PPG Beat
%   Detection (see README).

FS_MIN_DEFAULT = 30;                       % Hz; also the lowest rate allowed
if nargin < 3 || isempty(opts), opts = struct(); end
fs_min = FS_MIN_DEFAULT;
if isfield(opts, 'fs_min') && ~isempty(opts.fs_min)
    fs_min = opts.fs_min;
    if fs_min < FS_MIN_DEFAULT
        error('minel_ppg:fs_min', 'opts.fs_min must be at least %g Hz.', FS_MIN_DEFAULT);
    end
end
x = ppg(:);
x(~isfinite(x)) = 0;

if isfield(opts, 'stream') && strcmp(opts.stream, 'uniform')
    refine = isfield(opts, 'refine') && opts.refine;
    beats = detect_uniform(x, fs, refine);
    beats_int = beats;
    return;
end

% 1. decimation, leaving at least fs_min; optional block-average anti-aliasing
aa = 0;
if isfield(opts, 'aa') && ~isempty(opts.aa), aa = opts.aa; end
D = decimation_stride(fs, fs_min);
fsu = fs / D;
off = 0;                                       % filter delay in native samples
if D == 1 || aa == 0
    x = x(1:D:end);
else
    for p = 1:aa-1                             % extra boxcar stages at the native rate
        c = cumsum([0; x]);
        x = (c(D+1:end) - c(1:end-D)) / D;     % forward moving mean, length n-D+1
    end
    nb = floor(numel(x) / D);
    x = mean(reshape(x(1:nb*D), D, nb), 1)';   % block mean = boxcar + decimate
    off = aa * (D - 1) / 2;
end

% 2. derivative + band-pass (filter design cached per sampling rate)
[b, a] = bandpass_coeffs(fsu);
d = filtfilt(b, a, [0; diff(x)]);

% 3. sparse extremum coding
[sidx, sval] = sparse_extrema(d, fsu);
beats = zeros(0, 1); beats_int = zeros(0, 1);
if isempty(sidx), return; end

% 4-6. detection, blanking, quality gate (positions in the sparse stream)
p = canonical_pos(detect_beats(sval, 10), sidx);
p = blank_refractory(p, sidx, sval, fsu);
if isempty(p), return; end
p = quality_gate(p, sidx, sval, fsu);
if isempty(p), return; end
pk = sidx(p);

% 7. parabolic refinement of the derivative maximum, then native index
pkf = double(pk(:));
in = pkf > 1 & pkf < numel(d);                 % interior beats only
k  = pkf(in);
ym = d(k-1); y0 = d(k); yp = d(k+1);
den = ym - 2*y0 + yp;
delta = zeros(size(k));
ok = den < 0;                                  % true local maximum
delta(ok) = 0.5 * (ym(ok) - yp(ok)) ./ den(ok);
pkf(in) = k + max(-0.5, min(0.5, delta));

beats     = unique(round(D * (pkf - 1) + 1 + off));
beats_int = unique(round(D * (double(pk(:)) - 1) + 1 + off));
end


% =========================================================================
function [b, a] = bandpass_coeffs(fs)
% 0.5-5 Hz first-order Butterworth band-pass, cached for the last sampling rate.
persistent cached_fs cached_b cached_a
if isempty(cached_fs) || cached_fs ~= fs
    [cached_b, cached_a] = butter(1, [0.5 5] / (fs/2), 'bandpass');
    cached_fs = fs;
end
b = cached_b; a = cached_a;
end

% =========================================================================
function D = decimation_stride(fs, fs_min)
% Largest integer stride D that keeps fs/D >= fs_min (1 if none).
D = max(1, floor(fs / fs_min));
end

% =========================================================================
function p = canonical_pos(p, sidx)
% Sorted column of stream positions, each moved to the first occurrence of
% its sample index. sidx is sorted but can hold one index twice (a constant
% block has min and max at the same sample; a block has only two extrema,
% so never more than two copies). This reproduces ismember(sidx(p), sidx).
sidx = sidx(:);
p = sort(p(:));
if isempty(p), return; end
q = max(p - 1, 1);
dup = p > 1 & sidx(q) == sidx(p);
p(dup) = p(dup) - 1;
end

% =========================================================================
function [sidx, sval] = sparse_extrema(d, fs)
% Minimum and maximum of each 0.2 s block, with their sample indices.
% Block boundaries are round(0.2*j*fs), so non-integer block lengths
% (12.8 samples at 64 Hz, 6.25 at 31.25 Hz) introduce no drift.
T  = 10;                                   % values per second
sidx = []; sval = [];
if fs == round(fs)                         % integer rate: integer arithmetic
    FsI = fs;
    nb = floor(numel(d) * T / (2*FsI));    % number of blocks
    if nb < 20, return; end
    if mod(2*FsI, T) == 0                  % integer block length
        blk = 2*FsI / T;
        B = reshape(d(1:nb*blk), blk, []);
        [~, imin] = min(B, [], 1); [~, imax] = max(B, [], 1);
        sidx = sort([(0:nb-1)*blk + imin, (0:nb-1)*blk + imax]);
        sval = d(sidx);
        return;
    end
    e = floor((4*(0:nb)*FsI + T) / (2*T)); % = round(j*2*fs/T), exact
else                                       % non-integer rate (after decimation)
    nb = floor(numel(d) * T / (2*fs));
    if nb < 20, return; end
    e = round((0:nb) * (2*fs/T));
end
L = max(diff(e));
I = e(1:nb) + (1:L)';
I = min(I, e(2:nb+1));                     % clamp short blocks to their last sample
V = d(I);
[~, imin] = min(V, [], 1); [~, imax] = max(V, [], 1);
sidx = sort([e(1:nb) + imin, e(1:nb) + imax]);
sval = d(sidx);
end

% =========================================================================
function peaks = detect_uniform(x, fs, refine)
% Ablation: the same detector fed a uniformly resampled 10 Hz stream.
[b, a] = bandpass_coeffs(fs);
d = filtfilt(b, a, [0; diff(x)]);
[p, q] = rat(10 / fs);
y = resample(d, p, q);                     % anti-aliased, delay-compensated
peaks = zeros(0, 1);
if numel(y) < 40, return; end
sidx = round((0:numel(y)-1)' * fs / 10) + 1;
keep = sidx <= numel(d);
sidx = sidx(keep); sval = y(keep);
pos = canonical_pos(detect_beats(sval, 10), sidx);
pos = blank_refractory(pos, sidx, sval, fs);
pos = quality_gate(pos, sidx, sval, fs);
if isempty(pos), return; end
peaks = sidx(pos);
if refine                                  % +/-50 ms native-rate refinement
    w  = round(0.05 * fs);
    ii = min(max(peaks(:) + (-w:w), 1), numel(d));
    [~, j] = max(d(ii), [], 2);
    peaks = unique(ii(sub2ind(size(ii), (1:size(ii,1))', j)));
end
end

% =========================================================================
function keep = quality_gate(p, sidx, sval, fs)
% Rejects candidates whose upstroke slope is below 2.5% of four median
% absolute deviations of the 30 neighbouring candidates, and candidates
% below 2.5% of the 70th-percentile amplitude.
% p: candidate positions in the stream (sorted, canonical). Returns positions.
SQI_THRESH = 0.025;
keep = zeros(0, 1);
if isempty(p), return; end
sidx = sidx(:); sval = sval(:); pos = p(:);
n = numel(pos);
slopes = zeros(n, 1);
interior = pos > 1;
k = pos(interior);                         % upstroke from the previous extremum
slopes(interior) = max(0, (sval(k) - sval(k-1)) ./ ((sidx(k) - sidx(k-1)) / fs));
k = pos(~interior);                        % first stream value: use the next one
slopes(~interior) = abs(sval(k) - sval(k+1)) ./ ((sidx(k+1) - sidx(k)) / fs);
w = 30;
if n >= w
    m = movmad(slopes, w, 1); m(m < 1e-6) = 1e-6;
    sqi = min(1, slopes ./ (4*m));
else
    m = median(abs(slopes - median(slopes)));
    if m < 1e-6, m = 1e-6; end
    sqi = min(1, slopes / (4*m));
end
keep = pos(sqi >= SQI_THRESH);
if ~isempty(keep)
    h = sval(keep);
    keep = keep(h >= 0.025 * prctile(h, 70));
end
end

% =========================================================================
function keep = blank_refractory(p, sidx, sval, fs)
% Adaptive refractory blanking: hard refractory 0.27 s, adaptive refractory
% 0.45 median RR, dicrotic test (0.40-0.65 median RR, never earlier than
% 0.35 s, slope below 0.40 of the running median), and swap rescue of a
% candidate more than 1.5 times stronger than the beat blocking it.
% p: candidate positions in the stream (sorted, canonical). Returns positions.
SWAP_RATIO   = 1.5;
p = p(:);
n = numel(p);
if n == 0, keep = zeros(0, 1); return; end
cand = sidx(p); cand = cand(:);            % sample indices (timing)
amp  = sval(p); amp  = amp(:);
MIN_RR       = max(2, round(0.27 * fs));   % hard refractory (about 222 bpm)
MAX_REFR     = round(1.00 * fs);
ADAPT_FRAC   = 0.45;
MUST_FRAC    = 1.10;
PHASE_LO     = 0.40;
PHASE_HI     = 0.65;
AMP_RATIO    = 0.40;
AMP_MINB     = 5;
RR_WIN       = 8;
MIN_RR_FLOOR = round(0.35 * fs);           % dicrotic window floor (about 171 bpm)
rr_hist  = zeros(RR_WIN, 1); rr_count  = 0;
amp_hist = zeros(RR_WIN, 1); amp_hist(1) = amp(1); amp_count = 1;
acc = false(n, 1); acc(1) = true;
last_i = 1; prev_peak = NaN;
rr_dirty = true; amp_dirty = true; rr_med = 0; amp_med = 0;   % cached medians
for i = 2:n
    dt = cand(i) - cand(last_i);
    if rr_count == 0
        med_rr = 2 * MIN_RR;               % initial median RR, 0.54 s
    else
        if rr_dirty, rr_med = median(rr_hist(1:rr_count)); rr_dirty = false; end
        med_rr = rr_med;
    end
    % --- rejection tests ---
    rej = dt < MIN_RR;                                          % hard refractory
    if ~rej && dt <= MUST_FRAC * med_rr
        refr = min(MAX_REFR, max(MIN_RR, round(ADAPT_FRAC * med_rr)));
        if dt < refr
            rej = true;                                         % adaptive refractory
        elseif amp_count >= AMP_MINB
            if amp_dirty, amp_med = median(amp_hist(1:amp_count)); amp_dirty = false; end
            med_amp = amp_med;
            lo = max(PHASE_LO * med_rr, MIN_RR_FLOOR);
            hi = PHASE_HI * med_rr;
            if med_amp > 0 && dt >= lo && dt <= hi && amp(i) < AMP_RATIO * med_amp
                rej = true;                                     % dicrotic notch
            end
        end
    end
    % --- swap rescue ---
    if rej && amp(i) > SWAP_RATIO * amp(last_i)
        if isnan(prev_peak) || (cand(i) - prev_peak) >= MIN_RR
            acc(last_i) = false; acc(i) = true;
            if ~isnan(prev_peak) && rr_count > 0
                rr_hist(rr_count) = cand(i) - prev_peak; rr_dirty = true;
            end
            if amp_count > 0, amp_hist(amp_count) = amp(i); amp_dirty = true; end
            last_i = i;
            continue;
        end
    end
    if rej, continue; end
    % --- accept ---
    acc(i) = true;
    prev_peak = cand(last_i);
    if rr_count < RR_WIN
        rr_count = rr_count + 1; rr_hist(rr_count) = dt;
    else
        rr_hist(1:end-1) = rr_hist(2:end); rr_hist(end) = dt;
    end
    if amp_count < RR_WIN
        amp_count = amp_count + 1; amp_hist(amp_count) = amp(i);
    else
        amp_hist(1:end-1) = amp_hist(2:end); amp_hist(end) = amp(i);
    end
    rr_dirty = true; amp_dirty = true;
    last_i = i;
end
keep = p(acc);
end

% =========================================================================
function beats = detect_beats(s, fs)
% Rhythm-aware detector on the sparse stream s (fs = 10 values per second).
% Returns positions in s of accepted beats.
%   - Elgendi-style blocks of interest with two window sets (normal / fast),
%     selected by a heart-rate regime tracked with hysteresis (130/110 bpm).
%   - Dicrotic candidates resolved as a contest: fit to the expected rhythm
%     times saturated slope; timing weight shrinks at coarse RR resolution or
%     irregular rhythm. Window 0.30-0.70 median RR (0.35-0.70 fast), never
%     earlier than 0.35 s. A losing candidate is kept if its slope reaches
%     0.6 of the recent median; outside the fast regime this threshold relaxes
%     to 0.35 as the RR coefficient of variation rises from 0.10 to 0.25.
%   - Gated RR history with relock, searchback, final close-pair resolution.
beats = [];
if numel(s) < 15, return; end

AMP_CAP     = 1.5;        % slope advantage saturates here (x median slope)
AMP_W       = 0.3;        % slope weight when timing fully trusted
RR_GATE     = [0.6 1.6];  % plausible RR band for history update (x median)
GATE_ESCAPE = 3;          % consecutive out-of-band RRs -> relock
HR_INIT     = 120;        % initial regime threshold (bpm)
HR_UP       = 130;        % enter fast regime above
HR_DN       = 110;        % leave fast regime below

s = s(:)';
s_abs = abs(s);
N = numel(s);

% --- candidate masks (both window sets) and heart-rate regime ---
alpha = 0.202; beta = 0.0010;
off    = beta * mean(s_abs);
m_slow = (alpha*s_abs + (1-alpha)*conv(s_abs, ones(1,2)/2, 'same')) > ...
          conv(s_abs, ones(1,5)/5, 'same') + off;
m_fast = s_abs > conv(s_abs, ones(1,3)/3, 'same') + off;

fast = false(1, N);
st = find(diff([false, m_slow]) == 1);          % region starts = beat proxy
if numel(st) >= 4
    hr  = 60*fs ./ movmedian(diff(st), 9);
    reg = false(size(hr));
    state = hr(1) > HR_INIT;
    for k = 1:numel(hr)                         % hysteresis
        if state && hr(k) < HR_DN, state = false;
        elseif ~state && hr(k) > HR_UP, state = true; end
        reg(k) = state;
    end
    mk = false(1, N); mk(st(2:end)) = true;
    fast = reg(min(cumsum(mk) + 1, numel(reg)));
end

mask = m_slow; mask(fast) = m_fast(fast);
edges  = diff([false, mask, false]);
starts = find(edges == 1); ends = find(edges == -1) - 1;
if isempty(starts), return; end
cand_locs = zeros(size(starts));
cand_amps = zeros(size(starts));
for i = 1:numel(starts)
    [cand_amps(i), idx] = max(s(starts(i):ends(i)));
    cand_locs(i) = starts(i) + idx - 1;
end

% --- base amplitude threshold ---
win_len = min(N, round(58*fs));
base_thresh = 0.001 * prctile(s(end-win_len+1:end), 71);
valid = cand_amps > base_thresh;
cand_locs = cand_locs(valid);
cand_amps = cand_amps(valid);
if isempty(cand_locs), return; end
[cand_locs, ord] = sort(cand_locs);
cand_amps = cand_amps(ord);

% --- RR-adaptive acceptance with the dicrotic contest ---
min_dist    = 1;
rr_tol      = 0.123;
low_factor  = 0.366;
high_factor = 1.220;
MIN_RR_FLOOR = 0.35 * fs;                        % dicrotic window floor, 0.35 s (about 171 bpm)

n = numel(cand_locs);
cand_fast = fast(cand_locs);
qi = zeros(1, n);                                % accepted beat positions
qa = zeros(1, n);                                % accepted beat amplitudes
accepted = 0;
rr_hist  = zeros(1, 10); rr_cnt  = 0; rr_p  = 0;  % circular buffers
amp_hist = zeros(1, 10); amp_cnt = 0; amp_p = 0;
last_rr_added = false;
n_gate_miss   = 0;
dev = 0;                                          % running |RR - median|
% medians are recomputed only when their buffer changes (most candidates
% are rejected without changing either buffer)
rr_dirty = true; cv_dirty = true; amp_dirty = true;
rr_med = 0; cv_rr = 0; amp_med = 0;

for i = 1:n
    cur_t = cand_locs(i);
    cur_a = cand_amps(i);
    have_rr = accepted >= 2 && rr_cnt >= 1;

    if have_rr && rr_dirty
        rr_med = median(rr_hist(1:rr_cnt)); rr_dirty = false;
    end
    % irregularity-scaled threshold (not in the fast regime)
    if have_rr && rr_cnt >= 5 && ~cand_fast(i)
        if cv_dirty
            rr_mad = 1.4826 * median(abs(rr_hist(1:rr_cnt) - rr_med));
            cv_rr  = rr_mad / max(rr_med, eps); cv_dirty = false;
        end
        relax  = min(1, max(0, (cv_rr - 0.10) / 0.15));
    else
        relax  = 0;
    end
    notch_thr = 0.60 - 0.25 * relax;            % 0.60 (regular) -> 0.35 (irregular)

    if have_rr
        med_rr = rr_med;
        rr_dev = abs(cur_t - (qi(accepted) + med_rr)) / med_rr;
    else
        rr_dev = NaN;
    end
    mult = low_factor;
    if ~isnan(rr_dev) && rr_dev > rr_tol, mult = high_factor; end

    if accepted >= 1 && (cur_t - qi(accepted)) < min_dist, continue; end
    if cur_a <= mult*base_thresh, continue; end

    if have_rr
        dt = cur_t - qi(accepted);
        if cand_fast(i)
            lo = max(0.35 * med_rr, MIN_RR_FLOOR);
        else
            lo = max(0.30 * med_rr, MIN_RR_FLOOR);
        end
        hi = 0.70 * med_rr;
        if dt >= lo && dt <= hi                  % dicrotic contest
            if amp_dirty, amp_med = median(amp_hist(1:amp_cnt)); amp_dirty = false; end
            ref_a  = amp_med;
            prev_t = qi(accepted-1);
            cv  = 1.25*dev / med_rr;
            tw  = min(1, max(0, (med_rr - 4)/4)) * min(1, max(0, (0.25 - cv)/0.15));
            w   = AMP_W*tw + (1 - tw);
            cap = AMP_CAP*tw + 10*(1 - tw);
            sg  = max(1, 0.2*med_rr);
            dl  = (qi(accepted) - prev_t - med_rr) / sg;
            dc  = (cur_t - prev_t - med_rr) / sg;
            s_last = exp(-0.5*dl*dl)^tw * min(max(qa(accepted)/ref_a, 0.05), cap)^w;
            s_cur  = exp(-0.5*dc*dc)^tw * min(max(cur_a/ref_a, 0.05), cap)^w;
            if s_cur > s_last                    % current fits better: replace last
                qi(accepted) = cur_t;
                qa(accepted) = cur_a;
                if last_rr_added, rr_hist(rr_p) = cur_t - prev_t; rr_dirty = true; cv_dirty = true; end
                amp_hist(amp_p) = cur_a; amp_dirty = true;
                continue;
            end
            if cur_a / ref_a < notch_thr, continue; end   % weak loser: reject
        end
    end

    % --- accept ---
    accepted = accepted + 1;
    qi(accepted) = cur_t;
    qa(accepted) = cur_a;
    amp_p = mod(amp_p, 10) + 1; amp_hist(amp_p) = cur_a; amp_cnt = min(amp_cnt + 1, 10);
    amp_dirty = true;
    last_rr_added = false;
    if accepted >= 2
        new_rr = qi(accepted) - qi(accepted-1);
        ok = true;
        if rr_cnt >= 3
            ok = new_rr >= RR_GATE(1)*med_rr && new_rr <= RR_GATE(2)*med_rr;
            if ok
                n_gate_miss = 0;
            else
                n_gate_miss = n_gate_miss + 1;
                if n_gate_miss >= GATE_ESCAPE    % persistent change: relock
                    rr_cnt = 0; rr_p = 0; n_gate_miss = 0; ok = true; rr_dirty = true; cv_dirty = true;
                end
            end
        end
        if ok
            if rr_cnt == 0
                dev = 0.15*new_rr;
            else
                dev = 0.8*dev + 0.2*abs(new_rr - med_rr);
            end
            rr_p = mod(rr_p, 10) + 1; rr_hist(rr_p) = new_rr; rr_cnt = min(rr_cnt + 1, 10);
            last_rr_added = true; rr_dirty = true; cv_dirty = true;
        end
    end
end
qi = qi(1:accepted);
qa = qa(1:accepted);

% --- searchback: recover missed beats in long gaps ---
if accepted >= 2
    rr_all = diff(qi);
    med_rr = median(rr_all(max(1, end-9):end));
    ref_a  = median(qa);
    m_all  = median(rr_all);
    cv_g   = 1.4826*median(abs(rr_all - m_all)) / max(m_all, eps);
    idx = 2;
    while idx <= accepted
        a = qi(idx-1); b = qi(idx);
        mult = 1.5; if fast(a), mult = 1.8; end
        if b - a > mult*med_rr
            in_gap = cand_locs >= a+min_dist & cand_locs <= b-min_dist;
            cs = cand_locs(in_gap); ca = cand_amps(in_gap);
            if ~isempty(cs)
                tw  = min(1, max(0, (med_rr - 4)/4)) * min(1, max(0, (0.25 - cv_g)/0.15));
                w   = AMP_W*tw + (1 - tw);
                cap = AMP_CAP*tw + 10*(1 - tw);
                sg  = max(1, 0.2*med_rr);
                sc  = exp(-0.5*((cs - a - med_rr)/sg).^2).^tw .* ...
                      min(max(ca/ref_a, 0.05), cap).^w;
                [~, ic] = max(sc);
                if ca(ic) > 0.5*base_thresh
                    qi = [qi(1:idx-1), cs(ic), qi(idx:end)];
                    qa = [qa(1:idx-1), ca(ic), qa(idx:end)];
                    accepted = accepted + 1;
                    rr_all = diff(qi);
                    med_rr = median(rr_all(max(1, end-9):end));
                    continue;
                end
            end
        end
        idx = idx + 1;
    end
end

% --- final close-pair resolution ---
if accepted >= 3
    rr_all   = diff(qi);
    med_rr   = median(rr_all);
    refrac_f = max(round(0.25*fs), round(0.35*med_rr));
    refrac_s = max(round(0.30*fs), round(0.40*med_rr));
    ref_a    = median(qa);
    cv_g = 1.4826*median(abs(rr_all(:) - med_rr)) / max(med_rr, eps);
    tw  = min(1, max(0, (med_rr - 4)/4)) * min(1, max(0, (0.25 - cv_g)/0.15));
    w   = AMP_W*tw + (1 - tw);
    cap = AMP_CAP*tw + 10*(1 - tw);
    sg  = max(1, 0.2*med_rr);
    i = 2;
    while i <= numel(qi)
        if fast(qi(i)), refrac = refrac_f; else, refrac = refrac_s; end
        if (qi(i) - qi(i-1)) < refrac
            if i > 2
                anc = qi(i-2);
                d1 = (qi(i-1) - anc - med_rr) / sg;
                d2 = (qi(i)   - anc - med_rr) / sg;
                s1 = exp(-0.5*d1*d1)^tw * min(max(qa(i-1)/ref_a, 0.05), cap)^w;
                s2 = exp(-0.5*d2*d2)^tw * min(max(qa(i)/ref_a,   0.05), cap)^w;
                drop_second = s2 < s1;
            else
                drop_second = qa(i) < qa(i-1);
            end
            if drop_second
                qi(i) = []; qa(i) = [];
            else
                qi(i-1) = []; qa(i-1) = [];
                if i > 2, i = i - 1; end
            end
        else
            i = i + 1;
        end
    end
end
beats = qi;
end
