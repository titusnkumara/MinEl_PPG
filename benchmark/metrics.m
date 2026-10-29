function M = metrics()
%METRICS  Evaluation measures shared by run_paper.m and run_statistics.m.
%   M = metrics() returns handles to the functions below, e.g.
%   m = M.compute(ref, det, fs, tol, cpu, dur, cfg).
%
%   compute        per-recording sensitivity, PPV, F1, HR, SDNN, RMSSD, RR RMSE
%   aggregate      summary over recordings (medians, IQRs, coverage, timing)
%   print_table    console table: detectors x reference columns
%   clean_rr       artefact editing of RR intervals (0.25-2.0 s, 3 local MADs)
%   lag            per-detection pulse-transit offset (windowed median)
%   match_ref      logical: which reference beats have a detection within tol
%   match_mnn      mutual-nearest-neighbour matching within tol
%   interval_errors  RR-interval errors on consecutively matched beats (ms)
%   lag_window     offset window per database (60 s DaLiA, 300 s otherwise)
%   group_names    subgroup labels per database
M.compute         = @compute_metrics;
M.aggregate       = @aggregate_metrics;
M.print_table     = @print_table;
M.clean_rr        = @clean_rr_intervals;
M.lag             = @est_lag_per_detection;
M.match_ref       = @match_ref_logical;
M.match_mnn       = @match_mnn;
M.interval_errors = @interval_errors;
M.lag_window      = @lag_window;
M.group_names     = @group_names;
end


% =========================================================================
function m = compute_metrics(ref_peaks, det_peaks, Fs, mtol, cpu, dur, cfg)
m = struct('hr_ref',NaN,'hr_ppg',NaN,'sdnn_ref',NaN,'sdnn_ppg',NaN, ...
    'rmssd_ref',NaN,'rmssd_ppg',NaN,'rr_rmse',NaN, ...
    'sens',0,'ppv',0,'f1',0, ...
    'nref',numel(ref_peaks),'ndet',numel(det_peaks), ...
    'cpu',cpu,'dur',dur);
if isempty(ref_peaks) || isempty(det_peaks), return; end

% pulse-transit offset, removed equally for every detector
lag = est_lag_per_detection(ref_peaks, det_peaks, Fs, cfg);
det_peaks = round(det_peaks - lag);
det_peaks = det_peaks(det_peaks >= 1);
det_peaks = sort(unique(det_peaks));

[ref_m, ppg_m] = match_peaks(ref_peaks, det_peaks, Fs, mtol);
if numel(ref_m) < 3 || numel(ppg_m) < 3, return; end

rr_ref = diff(ref_m)/Fs; rr_ppg = diff(ppg_m)/Fs;
[rr_ref_c, keep_r] = clean_rr_intervals(rr_ref);
[rr_ppg_c, keep_p] = clean_rr_intervals(rr_ppg);

if numel(rr_ref_c) >= 5 && numel(rr_ppg_c) >= 5
    m.hr_ref = 60/mean(rr_ref_c); m.hr_ppg = 60/mean(rr_ppg_c);
    m.sdnn_ref = std(rr_ref_c); m.sdnn_ppg = std(rr_ppg_c);
    m.rmssd_ref = sqrt(mean(diff(rr_ref_c).^2));
    m.rmssd_ppg = sqrt(mean(diff(rr_ppg_c).^2));

    % RR RMSE from aligned 4 Hz tachograms
    t_ref = (ref_m(1:end-1)+ref_m(2:end))/2/Fs;
    t_ppg = (ppg_m(1:end-1)+ppg_m(2:end))/2/Fs;
    t_ref = t_ref(keep_r); t_ppg = t_ppg(keep_p);
    rr_ref_c = rr_ref_c(:); rr_ppg_c = rr_ppg_c(:);
    [t_ref, ia] = unique(t_ref); rr_ref_c = rr_ref_c(ia);
    [t_ppg, ia] = unique(t_ppg); rr_ppg_c = rr_ppg_c(ia);
    fs_i = 4;
    tg = min(t_ref(1),t_ppg(1)) : 1/fs_i : max(t_ref(end),t_ppg(end));
    rrg = interp1(t_ref, rr_ref_c, tg, 'linear','extrap');
    rpg = interp1(t_ppg, rr_ppg_c, tg, 'linear','extrap');
    [xc,lags] = xcorr(rrg-mean(rrg), rpg-mean(rpg), round(0.8*fs_i));
    [~,ix] = max(abs(xc)); lag = lags(ix)/fs_i;
    t_ppg2 = t_ppg + lag;
    tc = min(t_ref(1),t_ppg2(1)) : 1/fs_i : max(t_ref(end),t_ppg2(end));
    rrA = interp1(t_ref, rr_ref_c, tc, 'linear');
    rrB = interp1(t_ppg2, rr_ppg_c, tc, 'linear');
    v = ~isnan(rrA) & ~isnan(rrB);
    if sum(v) >= 5
        m.rr_rmse = sqrt(mean((rrA(v)-rrB(v)).^2));
    end
end

[m.sens, m.ppv, m.f1] = detection_metrics(ref_peaks, det_peaks, Fs, mtol);
end

% =========================================================================
function a = aggregate_metrics(M)
M = M(~cellfun(@isempty, M));
a = struct('nvalid',0,'nbeat',0, ...
    'hr_mae',NaN,'hr_mae_std',NaN,'hr_mae_med',NaN,'hr_mae_q1',NaN,'hr_mae_q3',NaN, ...
    'hr_bias',NaN,'hr_bias_med',NaN,'hr_std',NaN,'hr_corr',NaN, ...
    'sdnn_rmse',NaN,'sdnn_abs_med',NaN, ...
    'rmssd_rmse',NaN,'rmssd_abs_med',NaN, ...
    'rr_mean',NaN,'rr_med',NaN, ...
    'sens',NaN,'sens_med',NaN,'ppv',NaN,'ppv_med',NaN, ...
    'f1',NaN,'f1_med',NaN,'f1_q1',NaN,'f1_q3',NaN, ...
    'cov',NaN,'cov_nd',0,'cov_nr',0,'cpu',0,'dur',0,'rtf',NaN);
if isempty(M), return; end
M = [M{:}];

bv = [M.nref] > 0;
a.nbeat = sum(bv);
if a.nbeat >= 1
    sv = [M(bv).sens]; pv = [M(bv).ppv]; fv = [M(bv).f1];
    a.sens = mean(sv);  a.ppv = mean(pv);  a.f1 = mean(fv);
    a.sens_med = median(sv,'omitnan');
    a.ppv_med  = median(pv,'omitnan');
    fvv = fv(~isnan(fv));
    if ~isempty(fvv)
        a.f1_med = median(fvv);
        q = quantile(fvv,[0.25 0.75]);
        a.f1_q1 = q(1); a.f1_q3 = q(2);
    end
    a.cov_nd = sum([M(bv).ndet]); a.cov_nr = sum([M(bv).nref]);
    a.cov = 100*a.cov_nd/max(a.cov_nr,eps);
end

v = ~isnan([M.hr_ref]) & ~isnan([M.hr_ppg]);
a.nvalid = sum(v);
a.cpu = sum([M.cpu]); a.dur = sum([M.dur]);
a.rtf = a.dur/max(a.cpu,eps);
if a.nvalid < 2, return; end

e  = [M(v).hr_ref] - [M(v).hr_ppg];
ae = abs(e);
a.hr_mae     = mean(ae);
a.hr_mae_std = std(ae);
a.hr_mae_med = median(ae);
q = quantile(ae,[0.25 0.75]);
a.hr_mae_q1 = q(1); a.hr_mae_q3 = q(2);
a.hr_bias     = mean(e);
a.hr_bias_med = median(e);
a.hr_std      = std(e);
a.hr_corr     = corr([M(v).hr_ref]', [M(v).hr_ppg]');

se = [M(v).sdnn_ref]  - [M(v).sdnn_ppg];
re = [M(v).rmssd_ref] - [M(v).rmssd_ppg];
a.sdnn_rmse     = sqrt(mean(se.^2));
a.rmssd_rmse    = sqrt(mean(re.^2));
a.sdnn_abs_med  = median(abs(se),'omitnan');
a.rmssd_abs_med = median(abs(re),'omitnan');

rr = [M(v).rr_rmse];
a.rr_mean = mean(rr,'omitnan'); a.rr_med = median(rr,'omitnan');
end

% =========================================================================
function print_table(title_txt, refs, dets, A)
nr = numel(refs); nd = numel(dets);
W = 24;
fprintf('\n%s\n', title_txt);
hdr = sprintf('%-26s|', 'det \\ ref');
for r = 1:nr, hdr = [hdr sprintf(' %*s |', W, refs{r})]; end %#ok<AGROW>
fprintf('%s\n', hdr);
fprintf('%s\n', repmat('-', 1, numel(hdr)));
spec = {
    'F1 mean'                 @(a) sprintf('%.3f', a.f1)
    'F1 med (IQR)'            @(a) sprintf('%.3f (%.3f-%.3f)', a.f1_med, a.f1_q1, a.f1_q3)
    'Sens mean / med'         @(a) sprintf('%.3f / %.3f', a.sens, a.sens_med)
    'PPV mean / med'          @(a) sprintf('%.3f / %.3f', a.ppv, a.ppv_med)
    'HR MAE mean+-SD (bpm)'   @(a) sprintf('%.3f+-%.3f', a.hr_mae, a.hr_mae_std)
    'HR MAE med (IQR)'        @(a) sprintf('%.3f (%.3f-%.3f)', a.hr_mae_med, a.hr_mae_q1, a.hr_mae_q3)
    'HR bias mean / med'      @(a) sprintf('%.3f / %.3f', a.hr_bias, a.hr_bias_med)
    'HR LoA [mean-based]'     @(a) sprintf('[%.2f, %.2f]', a.hr_bias-1.96*a.hr_std, a.hr_bias+1.96*a.hr_std)
    'HR Corr'                 @(a) sprintf('%.3f', a.hr_corr)
    'SDNN RMSE / med|e| (s)'  @(a) sprintf('%.3f / %.3f', a.sdnn_rmse, a.sdnn_abs_med)
    'RMSSD RMSE / med|e| (s)' @(a) sprintf('%.3f / %.3f', a.rmssd_rmse, a.rmssd_abs_med)
    'RR RMSE mean / med (s)'  @(a) sprintf('%.3f / %.3f', a.rr_mean, a.rr_med)
    'Coverage'                @(a) sprintf('%.1f%% (%d/%d)', a.cov, a.cov_nd, a.cov_nr)
    'CPU total (s)'           @(a) sprintf('%.3f', a.cpu)
    'RTF (x1000)'             @(a) sprintf('%.1f', a.rtf/1000)
    'Valid (hr/beat)'         @(a) sprintf('%d/%d', a.nvalid, a.nbeat)
    };
for m = 1:size(spec,1)
    fprintf('%-26s\n', spec{m,1});
    for dt = 1:nd
        line = sprintf(' %-9s|', upper(dets{dt}));
        for r = 1:nr
            line = [line sprintf(' %*s |', W, spec{m,2}(A{dt,r}))]; %#ok<AGROW>
        end
        fprintf('%s\n', line);
    end
end
end

% =========================================================================
function [rr_clean, kept_idx] = clean_rr_intervals(rr)
if isempty(rr) || length(rr) < 3
    rr_clean = rr; kept_idx = 1:length(rr); return;
end
valid = (rr >= 0.25) & (rr <= 2.0);
rr = rr(valid);
if isempty(rr), rr_clean = []; kept_idx = []; return; end
window_size = min(10, length(rr));
kept = false(size(rr));
for i = 1:length(rr)
    start_idx = max(1, i - floor(window_size/2));
    end_idx = min(length(rr), i + floor(window_size/2));
    local_vals = rr(start_idx:end_idx);
    med = median(local_vals);
    mad = median(abs(local_vals - med));
    if mad == 0
        mad = std(local_vals) / 1.4826;
        if mad == 0, mad = 0.02; end
    end
    if abs(rr(i) - med) <= 3 * mad, kept(i) = true; end
end
rr_clean = rr(kept);
temp_idx = find(valid);
kept_idx = temp_idx(kept);
end

% =========================================================================
function [ref_matched, ppg_matched] = match_peaks(ref_peaks, ppg_peaks, Fs, tol)
if isempty(ppg_peaks) || isempty(ref_peaks)
    ref_matched = []; ppg_matched = []; return;
end
matched_ref_idx = zeros(length(ppg_peaks), 1);
for i = 1:length(ppg_peaks)
    [~, idx] = min(abs(ref_peaks - ppg_peaks(i)));
    if abs(ref_peaks(idx) - ppg_peaks(i)) / Fs <= tol
        matched_ref_idx(i) = idx;
    else
        matched_ref_idx(i) = NaN;
    end
end
unique_ref = unique(matched_ref_idx(~isnan(matched_ref_idx)));
keep_ppg = false(length(ppg_peaks), 1);
for j = 1:length(unique_ref)
    ref_idx = unique_ref(j);
    candidates = find(matched_ref_idx == ref_idx);
    if length(candidates) == 1
        keep_ppg(candidates) = true;
    else
        dists = abs(ppg_peaks(candidates) - ref_peaks(ref_idx));
        [~, best] = min(dists);
        keep_ppg(candidates(best)) = true;
    end
end
final_ppg = ppg_peaks(keep_ppg);
final_ref_idx = matched_ref_idx(keep_ppg);
final_ref = ref_peaks(final_ref_idx);
[final_ppg, sort_idx] = sort(final_ppg);
final_ref = final_ref(sort_idx);
ref_matched = final_ref(:);
ppg_matched = final_ppg(:);
end

% =========================================================================
function [sens, ppv, f1] = detection_metrics(ref_peaks, ppg_peaks, Fs, tol)
if isempty(ref_peaks) || isempty(ppg_peaks)
    sens = 0; ppv = 0; f1 = 0; return;
end
matched_ref_idx = nan(length(ppg_peaks),1);
for i = 1:length(ppg_peaks)
    [~, idx] = min(abs(ref_peaks - ppg_peaks(i)));
    if abs(ref_peaks(idx) - ppg_peaks(i)) / Fs <= tol
        matched_ref_idx(i) = idx;
    end
end
tp = length(unique(matched_ref_idx(~isnan(matched_ref_idx))));
fp = length(ppg_peaks) - tp;
fn = length(ref_peaks) - tp;
sens = tp / (tp + fn);
ppv = tp / (tp + fp);
f1 = 2 * sens * ppv / (sens + ppv);
if isnan(f1), f1 = 0; end
end

% =========================================================================
function mu = match_ref_logical(ref, det, Fs, tol)
mu = false(numel(ref),1);
if isempty(ref) || isempty(det), return; end
for i = 1:numel(ref)
    [~, j] = min(abs(det - ref(i)));
    if abs(det(j) - ref(i))/Fs <= tol, mu(i) = true; end
end
end

% =========================================================================
function [ri, dj] = match_mnn(r, d, tol)
ri = zeros(0,1); dj = zeros(0,1);
if isempty(r) || isempty(d), return; end
j  = knnsearch(d, r);
i2 = knnsearch(r, d);
idx = (1:numel(r))';
ok = i2(j) == idx & abs(d(j) - r) <= tol;
ri = idx(ok); dj = j(ok);
end

% =========================================================================
function [e_ms, key] = interval_errors(r, d, tol, Fs)
e_ms = zeros(0,1); key = zeros(0,1);
[ri, dj] = match_mnn(r, d, tol);
if numel(ri) < 2, return; end
ok = diff(ri) == 1 & diff(dj) == 1;
e  = (d(dj(2:end)) - d(dj(1:end-1))) - (r(ri(2:end)) - r(ri(1:end-1)));
e_ms = e(ok) / Fs * 1000;
k1 = ri(1:end-1); key = k1(ok);
end

% =========================================================================
function lag = est_lag(ref, det, max_samp)
lag = 0;
if numel(ref) < 3 || numel(det) < 3, return; end
j = knnsearch(det, ref);
d = det(j) - ref;
d = d(abs(d) <= max_samp);
if ~isempty(d), lag = median(d); end
end

% =========================================================================
function lag_per_det = est_lag_per_detection(ref, det, Fs, cfg)
% Windowed median nearest-neighbour offset, warm-started from the previous
% window. Returns one lag (samples) per detection. Records shorter than two
% windows use a single global offset.
ref = ref(:); det = det(:);
n_det = numel(det);
lag_per_det = zeros(n_det, 1);
if n_det == 0 || isempty(ref), return; end

t_det = (det - 1) / Fs;
t_ref = (ref - 1) / Fs;
total_sec = max([t_det; t_ref]);

if total_sec < 2 * cfg.LAG_WIN_SEC
    lag_per_det(:) = est_lag(ref, det, round(0.8*Fs));
    return;
end

prev_lag_samp = est_lag(ref, det, round(cfg.LAG_SEARCH_INITIAL_S * Fs));
win_sec = cfg.LAG_WIN_SEC;
win_starts = 0:win_sec:(total_sec - win_sec/2);
n_wins = numel(win_starts);

for w = 1:n_wins
    t0 = win_starts(w);
    t1 = t0 + win_sec;
    if w == n_wins, t1 = total_sec + 1; end
    in_win = t_det >= t0 & t_det < t1;
    ref_in = t_ref >= t0 & t_ref < t1;
    if sum(in_win) == 0, continue; end
    p_win = det(in_win);
    r_win = ref(ref_in);
    curr_lag_samp = prev_lag_samp;
    if numel(p_win) >= cfg.LAG_MIN_MATCHES && numel(r_win) >= cfg.LAG_MIN_MATCHES
        if w == 1
            ladder = cfg.LAG_SEARCH_INITIAL_S;
        else
            ladder = cfg.LAG_SEARCH_LADDER_S;
        end
        for search_s = ladder
            search_samp = round(search_s * Fs);
            p_adj = p_win - prev_lag_samp;
            j = knnsearch(p_adj, r_win);
            d = p_adj(j) - r_win;
            d = d(abs(d) <= search_samp);
            if numel(d) >= cfg.LAG_MIN_MATCHES
                med = median(d);
                mad = median(abs(d - med));
                if mad > 0
                    d = d(abs(d - med) <= 3*mad);
                end
                if numel(d) >= cfg.LAG_MIN_MATCHES
                    curr_lag_samp = median(d) + prev_lag_samp;
                    break;
                end
            end
        end
    end
    lag_per_det(in_win) = curr_lag_samp;
    prev_lag_samp = curr_lag_samp;
end
end

% =========================================================================
function w = lag_window(name)
% PPG-DaLiA: 60 s (independent ECG/PPG clocks); all others: 300 s.
if startsWith(name, 'dalia'), w = 60; else, w = 300; end
end

% =========================================================================
function gn = group_names(name)
switch name
    case 'mimic_af', gn = {'non-AF', 'AF'};
    case {'mimic_test', 'mimic_train'}, gn = {'adult', 'neonate'};
    otherwise, gn = {};
end
end
