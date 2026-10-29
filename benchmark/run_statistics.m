function run_statistics(ACC, cfg)
%RUN_STATISTICS  Every statistical analysis reported in the paper.
%   run_statistics(ACC, cfg) prints, in order:
%     A+B  paired Wilcoxon, proposed vs MSPTD_PB and vs MSPTDfast without
%          pre-filter (Holm per family)                       -> Table IV, S1
%     C    beat-level McNemar                                 -> Table S6
%     D    Bland-Altman agreement
%     E    AF screening (NEC) from each beat source           -> section 3.7
%     J    ECG reference quality (RR_CV)
%     PRIMARY  4-test Holm family on the held-out set         -> Table II
%     K    HRV non-inferiority on matched beats               -> Fig. 3
%     ABL  representation ablation                            -> Table V
%     G    RR-interval precision                              -> Table S5, Fig. 6a
%     PLAINREF vs MSPTD_PB
%     H    detection at tight tolerances                      -> Fig. 4a, Table S2
%     I    phase of false detections (PPG-DaLiA)
%     BEATCOUNT  reference-agreement filter                   -> Fig. 4b, Table S3
%     F    HRV by reference agreement
%     AF   exclusive beats in the AF database                 -> section 4.4
%     LoA  heart-rate limits of agreement, held-out set       -> section 3.1
%     END-TO-END heart rate and HRV                           -> Table III
%   Random seeds are fixed (rng 2024) so the output is reproducible.
blocks = {
    'A+B, C, D, E, J  paired tests, McNemar, Bland-Altman, AF screening, reference quality', @() inferential_stats(ACC, cfg)
    'PRIMARY  4-test Holm family',          @() primary_family_holm(ACC)
    'K  HRV non-inferiority',               @() hrv_noninferiority(ACC)
    'ABL  representation ablation',         @() ablation_stats(ACC)
    'G  RR-interval precision',             @() interval_precision(ACC, cfg)
    'PLAINREF vs MSPTD_PB',                 @() pair_compare(ACC, 'plainref', 'msptd_pb')
    'H  tight-tolerance detection',         @() tight_detection(ACC, cfg)
    'I  false-positive phase',              @() fp_phase(ACC, cfg)
    'BEATCOUNT  Tol = 5%',                  @() beatcount_agreement(ACC, 0.05)
    'BEATCOUNT  Tol = 10%',                 @() beatcount_agreement(ACC, 0.10)
    'BEATCOUNT  Tol = 15%',                 @() beatcount_agreement(ACC, 0.15)
    'BEATCOUNT  Tol = 20%',                 @() beatcount_agreement(ACC, 0.20)
    'F  HRV by reference agreement',        @() hrv_by_agreement(ACC)
    'AF  exclusive beats',                  @() af_exclusive(ACC, cfg)
    'LoA  heart-rate limits of agreement',  @() loa_check(ACC, cfg)
    'END-TO-END heart rate and HRV',        @() hrv_full(ACC)
    };
for b = 1:size(blocks, 1)
    fprintf('\n##### [%d/%d] %s #####\n', b, size(blocks, 1), blocks{b, 1});
    blocks{b, 2}();
end
end


% =========================================================================
% A+B, C, D, E, J
% =========================================================================
function inferential_stats(ACC, cfg)
M = metrics();
rng(2024);
ALPHA = 0.05;
du_all = cellfun(@(d) find(strcmpi(d,'user'),1), {ACC.dets}, 'UniformOutput', false);
ARMS = {'msptd_pb','msptd'};
METRICS = {'f1','sens','ppv','hr_mae','sdnn_err','rmssd_err','rr_rmse'};
METLAB  = {'F1','Sens','PPV','HR MAE (bpm)','SDNN err (s)','RMSSD err (s)','RR RMSE (s)'};

fprintf('\n=== A+B. Paired Wilcoxon (USER vs MSPTD variants), Holm alpha=%.2f ===\n', ALPHA);
fprintf('%-20s %-14s %-9s %-14s %4s %9s %9s %8s %8s %8s %6s %9s %s\n', ...
    'DB','arm','ref','metric','n','USER','MSPTD','med-d','CI-lo','CI-hi','r_rb','p','rej');
for ia = 1:numel(ARMS)
    arm = ARMS{ia};
    ds_all = cellfun(@(d) find(strcmpi(d, arm), 1), {ACC.dets}, 'UniformOutput', false);
    fprintf('\n--- USER vs %s ---\n', upper(arm));
    tbl = {}; allp = [];
    for k = 1:numel(ACC)
        du = du_all{k}; ds = ds_all{k};
        if isempty(du) || isempty(ds), continue; end
        for r = 1:numel(ACC(k).refs)
            for m = 1:numel(METRICS)
                [u, s] = paired_metric(ACC(k).R, r, du, ds, METRICS{m});
                if numel(u) < 5, continue; end
                p  = signrank(u, s);
                ci = boot_median_diff(u, s, 10000);
                tbl(end+1,:) = {ACC(k).name, ACC(k).refs{r}, METLAB{m}, numel(u), ...
                    median(u), median(s), median(u - s), ci(1), ci(2), rank_biserial(u, s), p}; %#ok<AGROW>
                allp(end+1) = p; %#ok<AGROW>
            end
        end
    end
    if isempty(allp), continue; end
    rej = holm(allp, ALPHA);
    for i = 1:size(tbl,1)
        fprintf('%-20s %-14s %-9s %-14s %4d %9.3f %9.3f %8.3f %8.3f %8.3f %6.2f %9.2e %s\n', ...
            tbl{i,1}, arm, tbl{i,2}, tbl{i,3}, tbl{i,4}, tbl{i,5}, tbl{i,6}, ...
            tbl{i,7}, tbl{i,8}, tbl{i,9}, tbl{i,10}, tbl{i,11}, repmat('*',1,rej(i)));
    end
    fprintf('  [%s: %d of %d tests rejected by Holm]\n', upper(arm), sum(rej), numel(allp));
end

fprintf('\n=== C. Beat-level McNemar (USER vs each MSPTD variant, ref col 1) ===\n');
for ia = 1:numel(ARMS)
    arm = ARMS{ia};
    ds_all = cellfun(@(d) find(strcmpi(d, arm), 1), {ACC.dets}, 'UniformOutput', false);
    fprintf('\n--- USER vs %s ---\n', upper(arm));
    for k = 1:numel(ACC)
        P = ACC(k).P; Fs = ACC(k).Fs; tol = ACC(k).tol;
        cfg.LAG_WIN_SEC = M.lag_window(ACC(k).name);
        du = du_all{k}; ds = ds_all{k};
        if isempty(du) || isempty(ds), continue; end
        b = 0; c = 0;
        for i = 1:numel(P)
            if isempty(P(i).refs), continue; end
            ref = P(i).refs{1};
            dU = round(P(i).dets{du} - M.lag(ref, P(i).dets{du}, Fs, cfg));
            dM = round(P(i).dets{ds} - M.lag(ref, P(i).dets{ds}, Fs, cfg));
            mu = M.match_ref(ref, dU, Fs, tol);
            mm = M.match_ref(ref, dM, Fs, tol);
            b = b + sum(mu & ~mm); c = c + sum(~mu & mm);
        end
        if b + c > 0
            chi2 = (abs(b-c)-1)^2/(b+c);
            fprintf('%-20s USER-only=%6d  %s-only=%6d  chi2=%.1f  p=%.2e\n', ...
                ACC(k).name, b, upper(arm), c, chi2, 1 - chi2cdf(chi2, 1));
        end
    end
end

fprintf('\n=== D. Bland-Altman USER (ref col 1) ===\n');
for k = 1:numel(ACC)
    R = ACC(k).R; du = du_all{k};
    if isempty(du), continue; end
    hr_r=[]; hr_p=[]; sd_r=[]; sd_p=[]; rm_r=[]; rm_p=[];
    for i = 1:numel(R)
        m = R(i).M{1,du};
        if isempty(m) || isnan(m.hr_ref) || isnan(m.hr_ppg), continue; end
        hr_r(end+1)=m.hr_ref;    hr_p(end+1)=m.hr_ppg;    %#ok<AGROW>
        sd_r(end+1)=m.sdnn_ref;  sd_p(end+1)=m.sdnn_ppg;  %#ok<AGROW>
        rm_r(end+1)=m.rmssd_ref; rm_p(end+1)=m.rmssd_ppg; %#ok<AGROW>
    end
    ba_summary(hr_r, hr_p, sprintf('%s HR (bpm)',  ACC(k).name));
    ba_summary(sd_r, sd_p, sprintf('%s SDNN (s)',  ACC(k).name));
    ba_summary(rm_r, rm_p, sprintf('%s RMSSD (s)', ACC(k).name));
end

fprintf('\n=== E. AF screening (NEC from USER beats) ===\n');
k = find(strcmp({ACC.name}, 'mimic_af'), 1);
if ~isempty(k)
    g = arrayfun(@(p) p.group, ACC(k).P);
    du = du_all{k};
    if any(g==1) && any(g==2) && ~isempty(du)
        nec = nan(numel(ACC(k).P), 1);
        for i = 1:numel(ACC(k).P)
            nec(i) = nec_from_peaks(ACC(k).P(i).dets{du}, ACC(k).Fs);
        end
        y = (g == 2); v = ~isnan(nec);
        if sum(v) >= 10
            [auc, auc_ci, thr, sens, spec] = nec_roc(nec(v), y(v));
            fprintf('%-11s NEC AUROC=%.3f [%.3f %.3f]  Youden thr=%.1f (sens=%.2f spec=%.2f)  LOSO-AUC=%.3f (n=%d)\n', ...
                ACC(k).name, auc, auc_ci(1), auc_ci(2), thr, sens, spec, loso_logistic(nec(v), y(v)), sum(v));
        end
    end

    % E2: the same screen on MinEl_PPG, MSPTD_PB and ECG (Pan-Tompkins) beats
    A = ACC(k); P = A.P; n = numel(P);
    y  = arrayfun(@(p) p.group, P) == 2;  y = y(:);
    dp = find(strcmpi(A.dets, 'msptd_pb'), 1);
    nec = nan(n, 3);
    for i = 1:n
        nec(i,1) = nec_from_peaks(P(i).dets{du}, A.Fs);
        nec(i,2) = nec_from_peaks(P(i).dets{dp}, A.Fs);
        nec(i,3) = nec_from_peaks(P(i).refs{1},  A.Fs);   % refs{1} = Pan-Tompkins
    end
    v = all(~isnan(nec), 2); nec = nec(v,:); y = y(v);
    names = {'MinEl_PPG', 'MSPTD_PB', 'ECG (PT)'};
    auc = zeros(1,3);
    for c = 1:3, [~,~,~,auc(c)] = perfcurve(y, nec(:,c), true); end
    rng(2024); B = 5000; iA = find(y); iN = find(~y); d = nan(B,3);
    pairs = [1 2; 1 3; 2 3];
    for b = 1:B                               % paired, class-stratified bootstrap
        s = [iA(randi(numel(iA), numel(iA), 1)); iN(randi(numel(iN), numel(iN), 1))];
        a = zeros(1,3);
        for c = 1:3, [~,~,~,a(c)] = perfcurve(y(s), nec(s,c), true); end
        d(b,:) = a(pairs(:,1)) - a(pairs(:,2));
    end
    fprintf('\n=== E2. AF screening by beat source (n = %d, AF = %d) ===\n', sum(v), sum(y));
    for c = 1:3, fprintf('%-10s AUROC = %.3f\n', names{c}, auc(c)); end
    for p = 1:3
        ci = prctile(d(:,p), [2.5 97.5]);
        fprintf('%-10s - %-10s  dAUROC = %+.3f  95%% CI [%+.3f, %+.3f]\n', ...
            names{pairs(p,1)}, names{pairs(p,2)}, auc(pairs(p,1)) - auc(pairs(p,2)), ci);
    end
end

reference_quality(ACC);
end

% =========================================================================
function [u, s] = paired_metric(R, r, du, ds, field)
u = nan(numel(R),1); s = nan(numel(R),1);
for k = 1:numel(R)
    mu = R(k).M{r,du}; ms = R(k).M{r,ds};
    if isempty(mu) || isempty(ms), continue; end
    switch field
        case 'hr_mae'
            if isnan(mu.hr_ref) || isnan(mu.hr_ppg) || isnan(ms.hr_ref) || isnan(ms.hr_ppg), continue; end
            u(k) = abs(mu.hr_ref - mu.hr_ppg); s(k) = abs(ms.hr_ref - ms.hr_ppg);
        case {'sdnn_err','rmssd_err'}
            f = strtok(field,'_');
            if isnan(mu.([f '_ref'])) || isnan(mu.([f '_ppg'])), continue; end
            if isnan(ms.([f '_ref'])) || isnan(ms.([f '_ppg'])), continue; end
            u(k) = abs(mu.([f '_ref']) - mu.([f '_ppg']));
            s(k) = abs(ms.([f '_ref']) - ms.([f '_ppg']));
        otherwise
            if ~isfield(mu,field) || ~isfield(ms,field), continue; end
            if isnan(mu.(field)) || isnan(ms.(field)), continue; end
            u(k) = mu.(field); s(k) = ms.(field);
    end
end
v = ~isnan(u) & ~isnan(s); u = u(v); s = s(v);
end

function rej = holm(p, alpha)
n = numel(p); [ps, ord] = sort(p); rej = false(1, n);
for j = 1:n
    if ps(j) <= alpha/(n-j+1), rej(ord(j)) = true; else, break; end
end
end

function r = rank_biserial(u, s)
d = u - s; d = d(d ~= 0); n = numel(d);
if n < 2, r = NaN; return; end
[~, ir] = sort(abs(d));
Tpos = 0;
for i = 1:n, if d(ir(i)) > 0, Tpos = Tpos + i; end, end
r = (4*Tpos)/(n*(n+1)) - 1;
end

function ci = boot_median_diff(u, s, nboot)
d = u - s; n = numel(d); b = zeros(nboot,1);
for i = 1:nboot, b(i) = median(d(randi(n,n,1))); end
ci = quantile(b, [0.025 0.975]);
end

function ba_summary(xr, xp, label)
d = xr - xp;
if numel(d) < 5, fprintf('%s: too few subjects\n', label); return; end
bias = mean(d); sd = std(d); n = numel(d);
nboot = 5000; B = zeros(nboot,1); L1 = zeros(nboot,1); L2 = zeros(nboot,1);
for i = 1:nboot
    dd = d(randi(n,n,1));
    B(i) = mean(dd);
    L1(i) = B(i) - 1.96*std(dd); L2(i) = B(i) + 1.96*std(dd);
end
fprintf('%s: bias=%.4f [%.4f,%.4f]  LoA=[%.4f, %.4f]  CI: L1 [%.4f,%.4f]  L2 [%.4f,%.4f]\n', ...
    label, bias, quantile(B,[.025 .975]), bias-1.96*sd, bias+1.96*sd, ...
    quantile(L1,[.025 .975]), quantile(L2,[.025 .975]));
end

function nec = nec_from_peaks(peaks, Fs)
% Number of non-empty cells in a 25 ms histogram of RR against successive RR
% differences, averaged over 64-beat windows (Lian et al 2011).
WINDOW_BEATS = 64; bin_width = 25;
if numel(peaks) < WINDOW_BEATS + 1, nec = NaN; return; end
peaks = sort(peaks(:));
rr_ms = diff(peaks)/Fs*1000;
n_windows = floor(numel(rr_ms)/WINDOW_BEATS);
edges_rr = 0:bin_width:3000; edges_drr = -2000:bin_width:2000;
vals = zeros(n_windows,1);
for w = 1:n_windows
    rr_win = rr_ms((w-1)*WINDOW_BEATS + 1 : w*WINDOW_BEATS);
    N = histcounts2(rr_win(2:end), diff(rr_win), edges_rr, edges_drr);
    vals(w) = sum(N(:) > 0);
end
nec = mean(vals);
end

function [auc, auc_ci, thr, sens, spec] = nec_roc(x, y)
[fpr, tpr, t, auc] = perfcurve(y, x, true);
nboot = 5000; a = zeros(nboot,1);
for i = 1:nboot
    idx = randi(numel(y), numel(y), 1);
    [~,~,~,a(i)] = perfcurve(y(idx), x(idx), true);
end
auc_ci = quantile(a, [.025 .975]);
[~, ib] = max(tpr - fpr);
thr = t(ib); sens = tpr(ib); spec = 1 - fpr(ib);
end

function auc = loso_logistic(x, y)
n = numel(y); pred = zeros(n,1);
for i = 1:n
    tr = true(n,1); tr(i) = false;
    mdl = fitglm(x(tr), y(tr), 'Distribution','binomial', 'LikelihoodPenalty','jeffreys-prior');
    pred(i) = predict(mdl, x(i));
end
[~,~,~,auc] = perfcurve(y, pred, true);
end

% =========================================================================
function reference_quality(ACC)
% J. RR coefficient of variation of each reference's beats that at least two
% of the three references agree on (lower = smoother intervals).
M = metrics();
TOL_CONS = round(0.15 * 125);
RR_MIN_1 = round(0.33 * 125);
RR_MIN_2 = round(0.25 * 125);
RR_MAX   = round(2.0  * 125);
REF_NAMES = {'PT','MinEl','TiTen'};
DB_NAMES  = {'mimic_test','mimic_train','bidmc','mimic_af'};
fprintf('\n=== J. ECG reference quality (RR_CV of per-reference consensus subset) ===\n');
fprintf('Lower RR_CV = smoother cardiac intervals = better timing accuracy.\n');
for dbk = 1:numel(DB_NAMES)
    k = find(strcmp({ACC.name}, DB_NAMES{dbk}), 1);
    if isempty(k), continue; end
    A = ACC(k);
    if ~A.is_ecg || numel(A.refs) < 3, continue; end
    nS = numel(A.P); nr = 3;
    g = arrayfun(@(p) p.group, A.P);
    rr_cv = nan(nS, nr);
    for i = 1:nS
        rr_min = RR_MIN_1; if g(i) == 2, rr_min = RR_MIN_2; end
        beats = cell(1, nr);
        for r = 1:nr, beats{r} = sort(double(A.P(i).refs{r}(:))); end
        all_b = sort(unique(cat(1, beats{:})));
        cons  = false(size(all_b));
        for j = 1:numel(all_b)
            votes = 0;
            for r = 1:nr
                if any(abs(beats{r} - all_b(j)) <= TOL_CONS), votes = votes+1; end
            end
            cons(j) = votes >= 2;
        end
        cons_beats = all_b(cons);
        for r = 1:nr
            b = beats{r};
            if numel(b) < 5, continue; end
            in_c = false(size(b));
            for j = 1:numel(b), in_c(j) = any(abs(cons_beats - b(j)) <= TOL_CONS); end
            cb = b(in_c);
            if numel(cb) >= 5
                rr = diff(cb);
                rr_v = rr(rr >= rr_min & rr <= RR_MAX);
                if numel(rr_v) >= 5, rr_cv(i,r) = std(rr_v) / mean(rr_v); end
            end
        end
    end
    gn = {'all'}; gm = {true(nS,1)};
    sub = M.group_names(DB_NAMES{dbk});
    if ~isempty(sub) && any(g==1) && any(g==2)
        gn = [gn, sub];
        gm = [gm, {(g==1)', (g==2)'}];
    end
    for gs = 1:numel(gn)
        idx = find(gm{gs});
        if numel(idx) < 5, continue; end
        fprintf('\n%s | %s (n=%d)\n', DB_NAMES{dbk}, gn{gs}, numel(idx));
        fprintf('  %-8s  %10s\n', 'Ref', 'RR_CV mean');
        for r = 1:nr, fprintf('  %-8s  %10.4f\n', REF_NAMES{r}, mean(rr_cv(idx,r), 'omitnan')); end
        pairs = [1 2; 1 3; 2 3];
        fprintf('  Wilcoxon RR_CV pairwise:\n');
        for p = 1:size(pairs,1)
            a = rr_cv(idx, pairs(p,1)); b = rr_cv(idx, pairs(p,2));
            v = ~isnan(a) & ~isnan(b);
            if sum(v) < 5, continue; end
            pv = signrank(a(v), b(v));
            sig = ''; if pv < 0.05, sig = '*'; end
            if pv < 0.001, sig = '***'; elseif pv < 0.01, sig = '**'; end
            fprintf('    %s vs %s: %.4f vs %.4f  p=%.4f %s\n', REF_NAMES{pairs(p,1)}, ...
                REF_NAMES{pairs(p,2)}, mean(a(v)), mean(b(v)), pv, sig);
        end
    end
end
end

% =========================================================================
% PRIMARY FAMILY and END-TO-END HRV
% =========================================================================
function primary_family_holm(ACC)
% Primary analysis (Methods 2.5): held-out set, PT reference, patient level,
% F1 plus end-to-end heart-rate, SDNN and RMSSD errors; Holm over 4 tests.
ALPHA = 0.05; DB = 'mimic_large_test';
k = find(strcmp({ACC.name}, DB), 1);
if isempty(k), fprintf('primary family: %s not run.\n', DB); return; end
du = find(strcmpi(ACC(k).dets,'user'),1); ds = find(strcmpi(ACC(k).dets,'msptd_pb'),1);
r  = find(strcmpi(ACC(k).refs,'pt'),1);
H = hrv_full_core(ACC, DB);
[f1u, f1s] = paired_metric(ACC(k).R, r, du, ds, 'f1');
[hu, hs] = hrv_pair(H, 1); [su, ss] = hrv_pair(H, 2); [mu, ms] = hrv_pair(H, 3);
T = { 'F1',                     f1u, f1s, '%.3f';
      'Heart-rate error (bpm)', hu,  hs,  '%.3f';
      'SDNN error (ms)',        su,  ss,  '%.1f';
      'RMSSD error (ms)',       mu,  ms,  '%.1f' };
n = size(T,1); p = nan(1,n); ci = nan(n,2);
rng(2024);                           % same bootstrap stream as hrv_full for rows 2-4
for i = [2 3 4 1]
    d = T{i,2} - T{i,3}; b = zeros(10000,1);
    for t = 1:10000, b(t) = median(d(randi(numel(d), numel(d), 1))); end
    ci(i,:) = quantile(b, [0.025 0.975]); p(i) = signrank(T{i,2}, T{i,3});
end
rej = holm(p, ALPHA);
fprintf('\n=== PRIMARY FAMILY: %d-test Holm (held-out, PT, %s n=%d) ===\n', n, H.unit, numel(f1u));
for i = 1:n
    u = T{i,2}; s = T{i,3}; f = T{i,4};
    fprintf(['%-24s n=%3d  U ' f ' (' f '-' f ')  PB ' f ' (' f '-' f ')  d %+.3f [%+.3f, %+.3f]  p=%.2e %s\n'], ...
        T{i,1}, numel(u), median(u), quantile(u,[.25 .75]), median(s), quantile(s,[.25 .75]), ...
        median(u-s), ci(i,:), p(i), repmat('*',1,rej(i)));
end
end

function hrv_full(ACC)
% End-to-end HR/HRV: each detector's full beat series against the
% reference's full series, extra and missed beats included (Table III).
rng(2024);
DBS = {'mimic_large_test','mimic_test','mimic_train','mimic_af','bidmc','capnobase','dalia'};
LAB = {'HR (bpm)', 'SDNN (ms)', 'RMSSD (ms)'};
for k = 1:numel(DBS)
    H = hrv_full_core(ACC, DBS{k}); if isempty(H), continue; end
    fprintf('\n=== %s (%s) | intervals edited: USER %.1f%%, PB %.1f%% ===\n', ...
        DBS{k}, H.unit, median(H.ED(:,1),'omitnan'), median(H.ED(:,2),'omitnan'));
    fprintf('%-11s %5s %9s %9s %9s %18s %10s  %s\n', 'measure', 'n', 'med|e| U', 'med|e| M', 'med(U-M)', '95% CI', 'p', 'verdict');
    for m = 1:3
        [u, s] = hrv_pair(H, m); if numel(u) < 5, continue; end
        dlt = u - s; b = zeros(10000,1);
        for t = 1:10000, b(t) = median(dlt(randi(numel(dlt), numel(dlt), 1))); end
        ci = quantile(b, [0.025 0.975]); p = signrank(u, s);
        if ci(2) < 0, vd = 'superior';
        elseif m > 1 && ci(2) < H.margin, vd = sprintf('non-inferior (margin %.1f ms)', H.margin);
        elseif ci(1) > 0, vd = 'worse';
        else, vd = 'no difference';
        end
        fprintf('%-11s %5d %9.3f %9.3f %+9.3f   [%+7.3f, %+7.3f] %10.2e  %s\n', ...
            LAB{m}, numel(u), median(u), median(s), median(dlt), ci(1), ci(2), p, vd);
    end
end
end

function H = hrv_full_core(ACC, dbname)
% End-to-end absolute errors per unit (patients if a database has repeated
% patients, else recordings). H.E: units x 3 (HR bpm, SDNN ms, RMSSD ms) x 2 (USER, PB).
H = []; SCL = [1 1000 1000];
a = find(strcmp({ACC.name}, dbname), 1); if isempty(a), return; end
A = ACC(a); Fs = A.Fs; n = numel(A.P);
dd = [find(strcmpi(A.dets,'user'),1), find(strcmpi(A.dets,'msptd_pb'),1)];
ri = 1; if A.is_ecg, ri = find(strcmpi(A.refs,'pt'),1); end
E = nan(n, 3, 2); ED = nan(n, 2);
for i = 1:n
    hr = hrv_of(A.P(i).refs{ri}, Fs); if isempty(hr), continue; end
    for d = 1:2
        hd = hrv_of(A.P(i).dets{dd(d)}, Fs); if isempty(hd), continue; end
        E(i,:,d) = abs(hd(1:3) - hr(1:3)) .* SCL; ED(i,d) = 100 * hd(4);
    end
end
pid = arrayfun(@(p) string(p.pid), A.P(:)); [~, ~, g] = unique(pid);
if max(g) < n
    E2 = nan(max(g), 3, 2); ED2 = nan(max(g), 2);
    for q = 1:max(g)
        E2(q,:,:) = mean(E(g==q,:,:), 1, 'omitnan'); ED2(q,:) = mean(ED(g==q,:), 1, 'omitnan');
    end
    E = E2; ED = ED2; H.unit = 'patients';
else
    H.unit = 'recordings';
end
H.E = E; H.ED = ED; H.margin = 2 * 1000 / A.Fs_ref;
end

function [u, s] = hrv_pair(H, m)
u = H.E(:,m,1); s = H.E(:,m,2); v = ~isnan(u) & ~isnan(s); u = u(v); s = s(v);
end

function out = hrv_of(beats, Fs)
% [HR, SDNN, RMSSD, fraction edited] from a full beat series; [] if too short.
M = metrics();
out = []; b = unique(round(beats(:))); if numel(b) < 7, return; end
rr = diff(b) / Fs; [rc, keep] = M.clean_rr(rr); rc = rc(:); keep = keep(:);
if numel(rc) < 5, return; end
ok = diff(keep) == 1; drr = diff(rc);          % only truly successive intervals
rmssd = NaN; if any(ok), rmssd = sqrt(mean(drr(ok).^2)); end
out = [60 / mean(rc), std(rc), rmssd, 1 - numel(rc) / numel(rr)];
end

% =========================================================================
% K. HRV NON-INFERIORITY ON MATCHED BEATS
% =========================================================================
function hrv_noninferiority(ACC)
% Margin = two sample periods of the reference clock. Non-inferior if the
% upper bound of the 95% bootstrap CI of the median paired difference in
% absolute error is below the margin; superior if below zero.
rng(2024);
NBOOT = 10000; ALPHA = 0.05; SCALE = 1000;
FIELDS = {'sdnn_err', 'rmssd_err', 'rr_rmse'};
LABS   = {'SDNN err (ms)', 'RMSSD err (ms)', 'RR RMSE (ms)'};
ARMS   = {'msptd_pb','msptd'};
fprintf('\n=== K. HRV Non-Inferiority (USER vs MSPTD variants; margin = 2 x sample period) ===\n');
fprintf('NI declared if upper bound of 95%% CI < delta_ms.\n');
fprintf('%-18s %-14s %-10s %-14s %5s %8s %8s %8s %9s %8s %8s  %s\n', ...
    'DB','arm','ref','metric','n','med(U)','med(M)','med(U-M)','CI-lo','CI-hi','delta','verdict');
for ia = 1:numel(ARMS)
    arm = ARMS{ia};
    fprintf('\n--- USER vs %s ---\n', upper(arm));
    for k = 1:numel(ACC)
        A = ACC(k);
        if startsWith(A.name, 'dalia_') || contains(A.name, '@'), continue; end
        delta = 2 * (1000 / A.Fs_ref);
        du = find(strcmpi(A.dets, 'user'), 1); ds = find(strcmpi(A.dets, arm), 1);
        if isempty(du) || isempty(ds), continue; end
        for r = 1:numel(A.refs)
            for f = 1:numel(FIELDS)
                [u, s] = paired_metric(A.R, r, du, ds, FIELDS{f});
                u = u*SCALE; s = s*SCALE; n = numel(u);
                if n < 5, continue; end
                d = u - s; bm = zeros(NBOOT, 1);
                for b = 1:NBOOT, bm(b) = median(d(randi(n, n, 1))); end
                lo = quantile(bm, ALPHA/2); hi = quantile(bm, 1 - ALPHA/2);
                if hi < 0, verdict = 'SUPERIOR';
                elseif hi < delta, verdict = 'NON-INFERIOR';
                else, verdict = 'inconclusive';
                end
                fprintf('%-18s %-14s %-10s %-14s %5d %8.2f %8.2f %8.2f %9.2f %8.2f %8.2f  %s\n', ...
                    A.name, arm, A.refs{r}, LABS{f}, n, median(u), median(s), median(d), lo, hi, delta, verdict);
            end
        end
    end
end
end

% =========================================================================
% ABLATION and PAIRWISE ARM COMPARISON
% =========================================================================
function ablation_stats(ACC)
rng(2024); ALPHA = 0.05;
METRICS = {'f1','sens','ppv','hr_mae','rr_rmse'};
METLAB  = {'F1','Sens','PPV','HR MAE (bpm)','RR RMSE (s)'};
PAIRS = {'user','plain'; 'user','plainref'; 'plainref','plain'; 'plainref','msptd_pb'; 'user','msptd_pb'};
tbl = {}; allp = [];
for k = 1:numel(ACC)
    if startsWith(ACC(k).name, 'dalia_'), continue; end
    R = ACC(k).R; dets = ACC(k).dets;
    for q = 1:size(PAIRS,1)
        ia = find(strcmpi(dets, PAIRS{q,1}), 1); ib = find(strcmpi(dets, PAIRS{q,2}), 1);
        if isempty(ia) || isempty(ib), continue; end
        for m = 1:numel(METRICS)
            [u, s] = paired_metric(R, 1, ia, ib, METRICS{m});
            if numel(u) < 5, continue; end
            p = signrank(u, s);
            tbl(end+1,:) = {ACC(k).name, sprintf('%s-%s', PAIRS{q,1}, PAIRS{q,2}), METLAB{m}, ...
                numel(u), median(u), median(s), median(u - s), p}; %#ok<AGROW>
            allp(end+1) = p; %#ok<AGROW>
        end
    end
end
rej = holm(allp, ALPHA);
fprintf('\n=== ABLATION. Paired Wilcoxon between representations, Holm alpha=%.2f (%d tests) ===\n', ALPHA, numel(allp));
fprintf('%-12s %-16s %-14s %4s %9s %9s %9s %9s %s\n', 'DB','A - B','metric','n','med A','med B','med A-B','p','rej');
for i = 1:size(tbl,1)
    fprintf('%-12s %-16s %-14s %4d %9.3f %9.3f %9.3f %9.2e %s\n', ...
        tbl{i,1}, tbl{i,2}, tbl{i,3}, tbl{i,4}, tbl{i,5}, tbl{i,6}, tbl{i,7}, tbl{i,8}, repmat('*',1,rej(i)));
end
end

function pair_compare(ACC, armA, armB)
M = metrics();
rng(2024); ALPHA = 0.05; NBOOT = 10000;
METRICS = {'f1','sens','ppv','hr_mae','sdnn_err','rmssd_err','rr_rmse'};
METLAB  = {'F1','Sens','PPV','HR MAE (bpm)','SDNN err (s)','RMSSD err (s)','RR RMSE (s)'};
fprintf('\n=== Paired Wilcoxon: %s vs %s (Holm alpha=%.2f) ===\n', upper(armA), upper(armB), ALPHA);
fprintf('%-20s %-9s %-14s %4s %9s %9s %9s %9s %9s %6s %9s %s\n', ...
    'DB','ref','metric','n','A','B','med(A-B)','CI-lo','CI-hi','r_rb','p','rej');
tbl = {}; allp = [];
for k = 1:numel(ACC)
    A = ACC(k);
    if startsWith(A.name, 'dalia_') || contains(A.name, '@'), continue; end
    ia = find(strcmpi(A.dets, armA), 1); ib = find(strcmpi(A.dets, armB), 1);
    if isempty(ia) || isempty(ib), continue; end
    gn = M.group_names(A.name); g = [A.P.group]';
    labels = {'all'}; masks = {true(numel(A.R),1)};
    if numel(A.R) ~= numel(A.P), gn = {}; end
    for q = 1:numel(gn)
        if any(g == q), labels{end+1} = gn{q}; masks{end+1} = (g == q); end %#ok<AGROW>
    end
    for gs = 1:numel(labels)
        idx = find(masks{gs});
        if numel(idx) < 5, continue; end
        for r = 1:numel(A.refs)
            for m = 1:numel(METRICS)
                [u, s] = paired_metric(A.R(idx), r, ia, ib, METRICS{m});
                if numel(u) < 5, continue; end
                p  = signrank(u, s);
                ci = boot_median_diff(u, s, NBOOT);
                lab = A.name; if gs > 1, lab = sprintf('%s@%s', A.name, labels{gs}); end
                tbl(end+1,:) = {lab, A.refs{r}, METLAB{m}, numel(u), median(u), median(s), ...
                    median(u - s), ci(1), ci(2), rank_biserial(u, s), p}; %#ok<AGROW>
                allp(end+1) = p; %#ok<AGROW>
            end
        end
    end
end
if isempty(allp), fprintf('  (no valid pairs)\n'); return; end
rej = holm(allp, ALPHA);
for i = 1:size(tbl,1)
    fprintf('%-20s %-9s %-14s %4d %9.3f %9.3f %9.3f %9.3f %9.3f %6.2f %9.2e %s\n', ...
        tbl{i,1}, tbl{i,2}, tbl{i,3}, tbl{i,4}, tbl{i,5}, tbl{i,6}, ...
        tbl{i,7}, tbl{i,8}, tbl{i,9}, tbl{i,10}, tbl{i,11}, repmat('*',1,rej(i)));
end
end

% =========================================================================
% G. RR-INTERVAL PRECISION
% =========================================================================
function interval_precision(ACC, cfg)
% Offset-free RR-interval errors pooled over recordings, against Pan-Tompkins
% for ECG-referenced databases. Common sets compare arms on identical intervals.
M = metrics();
TOL = 0.1;
SETS = {{'user','plain','plainref'}, {'user','msptd_pb'}, {'user','msptd'}};
fprintf('\n=== G. RR-interval precision vs reference (ms; offset removed per recording) ===\n');
for k = 1:numel(ACC)
    A = ACC(k);
    if startsWith(A.name, 'dalia_') || contains(A.name, '@'), continue; end
    cfg.LAG_WIN_SEC = M.lag_window(A.name);
    Fs = A.Fs; nd = numel(A.dets); nS = numel(A.P);
    ri = 1;
    if A.is_ecg, ip = find(strcmpi(A.refs, 'pt'), 1); if ~isempty(ip), ri = ip; end, end
    es = cell(nS, nd); ks = cell(nS, nd); nRIs = zeros(nS, 1);
    for i = 1:nS
        r = sort(double(A.P(i).refs{ri}(:)));
        if numel(r) < 20, continue; end
        nRIs(i) = numel(r) - 1;
        for d = 1:nd
            dd = sort(double(A.P(i).dets{d}(:)));
            if numel(dd) < 20, continue; end
            dd = dd - M.lag(r, dd, Fs, cfg);
            [es{i,d}, ks{i,d}] = M.interval_errors(r, dd, round(TOL*Fs), Fs);
        end
    end
    gn = M.group_names(A.name); g = [A.P.group]';
    labels = {'all'}; masks = {true(nS, 1)};
    for q = 1:numel(gn), labels{end+1} = gn{q}; masks{end+1} = (g == q); end %#ok<AGROW>
    for q = 1:numel(labels)
        idx = find(masks{q});
        if numel(idx) < 3, continue; end
        fprintf('\n%s | %s (n = %d) | reference: %s | sample resolution %.1f ms\n', ...
            A.name, labels{q}, numel(idx), A.refs{ri}, 1000/Fs);
        fprintf('  %-12s %8s %9s %8s %8s %9s %9s\n', 'detector', 'n', 'coverage', '<=20ms', '<=50ms', 'med|e|', 'P95|e|');
        for d = 1:nd
            e = vertcat(es{idx,d});
            if isempty(e), continue; end
            fprintf('  %-12s %8d %8.1f%% %7.1f%% %7.1f%% %8.1f  %8.1f\n', upper(A.dets{d}), numel(e), ...
                100*numel(e)/max(1, sum(nRIs(idx))), 100*mean(abs(e) <= 20), 100*mean(abs(e) <= 50), ...
                median(abs(e)), prctile(abs(e), 95));
        end
        for sIdx = 1:numel(SETS)
            names = SETS{sIdx};
            f = cellfun(@(nm) find(strcmpi(A.dets, nm), 1), names, 'UniformOutput', false);
            if any(cellfun(@isempty, f)), continue; end
            f = cell2mat(f); nf = numel(f); C = cell(1, nf);
            for i = idx'
                if any(cellfun(@isempty, ks(i, f))), continue; end
                c = ks{i, f(1)};
                for q2 = 2:nf, c = intersect(c, ks{i, f(q2)}); end
                for q2 = 1:nf
                    [~, ia] = ismember(c, ks{i, f(q2)});
                    C{q2} = [C{q2}; es{i, f(q2)}(ia)];
                end
            end
            if ~isempty(C{1})
                fprintf('  common set (%s), n = %d:\n', strjoin(upper(names), ' / '), numel(C{1}));
                for q2 = 1:nf
                    fprintf('  %-12s %8s %9s %7.1f%% %7.1f%% %8.1f\n', upper(A.dets{f(q2)}), '', '', ...
                        100*mean(abs(C{q2}) <= 20), 100*mean(abs(C{q2}) <= 50), median(abs(C{q2})));
                end
            end
        end
    end
end
end

% =========================================================================
% H. DETECTION AT TIGHT TOLERANCES
% =========================================================================
function tight_detection(ACC, cfg)
% Offset removed per recording, mutual-nearest matching at 75 and 100 ms
% (the main tables use 150 ms). Held-out tests are at patient level.
M = metrics();
TOLS = [0.075 0.1];
fprintf('\n=== H. Detection at tight tolerance (offset removed per recording; mutual-nearest matching) ===\n');
stat = {};
for k = 1:numel(ACC)
    A = ACC(k);
    if startsWith(A.name, 'dalia_') || contains(A.name, '@'), continue; end
    cfg.LAG_WIN_SEC = M.lag_window(A.name);
    Fs = A.Fs; nd = numel(A.dets); nS = numel(A.P);
    ri = 1;
    if A.is_ecg, ip = find(strcmpi(A.refs, 'pt'), 1); if ~isempty(ip), ri = ip; end, end
    F = nan(nS, nd, numel(TOLS)); SE = F; PP = F;
    for i = 1:nS
        r = sort(double(A.P(i).refs{ri}(:)));
        if numel(r) < 20, continue; end
        for d = 1:nd
            dd = sort(double(A.P(i).dets{d}(:)));
            if numel(dd) < 20, continue; end
            dd = dd - M.lag(r, dd, Fs, cfg);
            for t = 1:numel(TOLS)
                [rm, ~] = M.match_mnn(r, dd, round(TOLS(t)*Fs));
                tp = numel(rm);
                SE(i,d,t) = tp / numel(r); PP(i,d,t) = tp / numel(dd);
                F(i,d,t) = 2*tp / (numel(r) + numel(dd));
            end
        end
    end
    gn = M.group_names(A.name); g = [A.P.group]';
    labels = {'all'}; masks = {true(nS, 1)};
    for q = 1:numel(gn), labels{end+1} = gn{q}; masks{end+1} = (g == q); end %#ok<AGROW>
    du = find(strcmpi(A.dets, 'user'), 1);
    pid = [A.P.pid]';
    for q = 1:numel(labels)
        idx = find(masks{q});
        if numel(idx) < 3, continue; end
        for t = 1:numel(TOLS)
            fprintf('\n%s | %s (n = %d) | reference: %s | tolerance %.2f s\n', ...
                A.name, labels{q}, numel(idx), A.refs{ri}, TOLS(t));
            fprintf('  %-12s %6s %6s %6s\n', 'detector', 'Sens', 'PPV', 'F1');
            for d = 1:nd
                fprintf('  %-12s %6.3f %6.3f %6.3f\n', upper(A.dets{d}), median(SE(idx,d,t), 'omitnan'), ...
                    median(PP(idx,d,t), 'omitnan'), median(F(idx,d,t), 'omitnan'));
                if ~isempty(du) && d ~= du
                    [~,~,gi] = unique(pid(idx));
                    u = accumarray(gi, F(idx,du,t), [], @(x) mean(x,'omitnan'));
                    s = accumarray(gi, F(idx,d,t),  [], @(x) mean(x,'omitnan'));
                    v = ~isnan(u) & ~isnan(s);
                    if sum(v) >= 5
                        stat(end+1,:) = {A.name, labels{q}, TOLS(t), sprintf('user vs %s', A.dets{d}), ...
                            sum(v), median(u(v)), median(s(v)), signrank(u(v), s(v))}; %#ok<AGROW>
                    end
                end
            end
        end
    end
end
if ~isempty(stat)
    p = cell2mat(stat(:,8)); rej = holm(p', 0.05);
    fprintf('\n--- paired Wilcoxon on per-recording F1 (USER vs other arm), Holm over %d tests ---\n', numel(p));
    for i = 1:size(stat, 1)
        fprintf('%-12s %-8s tol %.2f  %-22s n=%3d  %.3f vs %.3f  p=%.2e %s\n', stat{i,1}, stat{i,2}, ...
            stat{i,3}, stat{i,4}, stat{i,5}, stat{i,6}, stat{i,7}, stat{i,8}, repmat('*', 1, rej(i)));
    end
end
end

% =========================================================================
% I. PHASE OF FALSE DETECTIONS (PPG-DaLiA)
% =========================================================================
function fp_phase(ACC, cfg)
% Phase 0 = a reference beat, 1 = the next one. False detections near 0 are
% double detections; mid-cycle ones are diastolic or dicrotic events.
M = metrics();
names = {'dalia_sitting', 'dalia_working'};
TOL = 0.15; edges = 0:0.1:1;
fprintf('\n=== I. Phase of false detections within the reference cycle (DaLiA) ===\n');
for n = 1:numel(names)
    k = find(strcmp({ACC.name}, names{n}), 1);
    if isempty(k), continue; end
    A = ACC(k); Fs = A.Fs; nd = numel(A.dets); nS = numel(A.P);
    cfg.LAG_WIN_SEC = M.lag_window(A.name);
    fprintf('\n%s\n  %-12s %7s   share of false detections by phase bin (bin starts: %s)\n', ...
        names{n}, 'detector', 'FP/100', sprintf('%4.1f', edges(1:end-1)));
    for d = 1:nd
        ph = []; nfp = 0; nref = 0;
        for i = 1:nS
            r = unique(double(A.P(i).refs{1}(:)));
            dd = sort(double(A.P(i).dets{d}(:)));
            if numel(r) < 20 || numel(dd) < 20, continue; end
            nref = nref + numel(r);
            dd = dd - M.lag(r, dd, Fs, cfg);
            [~, dj] = M.match_mnn(r, dd, round(TOL*Fs));
            fp = true(numel(dd), 1); fp(dj) = false;
            t = dd(fp); nfp = nfp + numel(t);
            kk = discretize(t, r);
            ok = ~isnan(kk) & kk < numel(r);
            kk = kk(ok); t = t(ok);
            ph = [ph; (t - r(kk)) ./ (r(kk+1) - r(kk))]; %#ok<AGROW>
        end
        h = histcounts(ph, edges);
        if sum(h) == 0, continue; end
        fprintf('  %-12s %7.1f   ', upper(A.dets{d}), 100*nfp/max(1, nref));
        fprintf('%5.0f', 100*h/sum(h));
        fprintf('\n');
    end
end
end

% =========================================================================
% BEAT-COUNT AGREEMENT FILTER
% =========================================================================
function beatcount_agreement(ACC, Tol)
% Keep recordings whose three ECG references agree on beat count within a
% relative spread Tol: (max - min) / max <= Tol.
fprintf('\n=== BEAT-COUNT AGREEMENT FILTER (Tol = %.0f%%) ===\n', Tol*100);
fprintf('%-22s %6s %6s %6s %8s %8s %8s %10s %8s\n', ...
    'Database','n_full','n_kept','%','F1_USER','F1_MSPTD','dF1','p','status');
for k = 1:numel(ACC)
    A = ACC(k);
    if ~A.is_ecg || numel(A.refs) < 3, continue; end
    if startsWith(A.name,'dalia_') || contains(A.name,'@'), continue; end
    if numel(A.R) ~= numel(A.P)
        fprintf('%-22s  [skipped: patient-collapsed, cannot filter record-wise]\n', A.name);
        continue;
    end
    ip = find(strcmpi(A.refs,'pt'),1);
    du = find(strcmpi(A.dets,'user'),1);
    dm = find(strcmpi(A.dets,'msptd_pb'),1);
    if any(cellfun(@isempty,{ip,du,dm})), continue; end
    nS = numel(A.P);
    keep = false(nS,1);
    for i = 1:nS
        c = cellfun(@numel, A.P(i).refs(1:numel(A.refs)));
        keep(i) = numel(c) >= 3 && (max(c) - min(c)) / max(c) <= Tol;
    end
    fU = nan(nS,1); fM = nan(nS,1);
    for i = 1:nS
        MU = A.R(i).M{ip,du}; MM = A.R(i).M{ip,dm};
        if ~isempty(MU), fU(i) = MU.f1; end
        if ~isempty(MM), fM(i) = MM.f1; end
    end
    mk = keep & ~isnan(fU) & ~isnan(fM);
    nk = sum(mk);
    if nk < 5
        fprintf('%-22s %6d %6d %6s %8s %8s %8s %10s  too few\n', A.name, nS, nk, '-', '-', '-', '-', '-');
        continue;
    end
    fprintf('%-22s %6d %6d %5.0f%% %8.3f %8.3f %+8.3f %10.2e  ok\n', A.name, nS, nk, 100*nk/nS, ...
        median(fU(mk)), median(fM(mk)), median(fU(mk) - fM(mk)), signrank(fU(mk), fM(mk)));
end
end

% =========================================================================
% F. HRV BY REFERENCE AGREEMENT
% =========================================================================
function hrv_by_agreement(ACC)
fprintf('\n=== F. HRV stratified by ECG-reference agreement ===\n');
for k = 1:numel(ACC)
    if ~ACC(k).is_ecg || numel(ACC(k).refs) < 2, continue; end
    P = ACC(k).P; R = ACC(k).R;
    if numel(R) ~= numel(P), continue; end
    nr = numel(ACC(k).refs);
    du = find(strcmpi(ACC(k).dets,'user'),1); ds = find(strcmpi(ACC(k).dets,'msptd_pb'),1);
    if isempty(du) || isempty(ds), continue; end
    md = zeros(numel(P),1);
    for i = 1:numel(P)
        c = cellfun(@numel, P(i).refs(1:nr));
        md(i) = max([abs(c(1)-c(2)), abs(c(1)-c(3)), abs(c(2)-c(3))]) / max(c);
    end
    for thr = [0.05 0.10 1.0]
        idx = find(md <= thr);
        if numel(idx) < 5, continue; end
        e = nan(numel(idx), 4);                  % SDNN U, SDNN M, RMSSD U, RMSSD M
        for j = 1:numel(idx)
            Mu = R(idx(j)).M{1, du}; Ms = R(idx(j)).M{1, ds};
            if isempty(Mu) || isempty(Ms), continue; end
            if any(isnan([Mu.sdnn_ref Mu.sdnn_ppg Ms.sdnn_ref Ms.sdnn_ppg])), continue; end
            e(j,:) = abs([Mu.sdnn_ref - Mu.sdnn_ppg, Ms.sdnn_ref - Ms.sdnn_ppg, ...
                          Mu.rmssd_ref - Mu.rmssd_ppg, Ms.rmssd_ref - Ms.rmssd_ppg]);
        end
        e = e(~isnan(e(:,1)) & ~isnan(e(:,2)), :);
        if size(e,1) < 5, continue; end
        fprintf('%-11s ref-ok<=%.2f: n=%3d  SDNN %.3f/%.3f (p=%.3f)  RMSSD %.3f/%.3f (p=%.3f)\n', ...
            ACC(k).name, thr, numel(idx), mean(e(:,1)), mean(e(:,2)), signrank(e(:,1), e(:,2)), ...
            mean(e(:,3)), mean(e(:,4)), signrank(e(:,3), e(:,4)));
    end
end
end

% =========================================================================
% AF. EXCLUSIVE BEATS IN THE AF DATABASE
% =========================================================================
function af_exclusive(ACC, cfg)
% Beats found by only one detector, by preceding-RR ratio and by reference
% agreement, and sensitivity on beats all three references agree on.
M = metrics();
k = find(strcmp({ACC.name}, 'mimic_af'), 1);
if isempty(k), fprintf('mimic_af not run, skipped.\n'); return; end
A = ACC(k); Fs = A.Fs; tol = A.tol;
du = find(strcmpi(A.dets, 'user'), 1); dm = find(strcmpi(A.dets, 'msptd_pb'), 1);
cfg.LAG_WIN_SEC = M.lag_window(A.name);
edges = [0 0.6 0.8 1.0 1.2 1.5 Inf]; nb = numel(edges) - 1;
gn = M.group_names('mimic_af');
for grp = 1:2
    cntU = zeros(1, nb); cntM = zeros(1, nb); tot = zeros(1, nb);
    agU = [0 0]; agM = [0 0]; sU = []; sM = [];
    for i = find([A.P.group] == grp)
        r = A.P(i).refs{1}(:); if numel(r) < 10, continue; end
        dU = round(A.P(i).dets{du} - M.lag(r, A.P(i).dets{du}, Fs, cfg));
        dM = round(A.P(i).dets{dm} - M.lag(r, A.P(i).dets{dm}, Fs, cfg));
        mu = M.match_ref(r, dU, Fs, tol); mm = M.match_ref(r, dM, Fs, tol);
        rrp = [NaN; diff(r)]; ratio = rrp ./ movmedian(rrp, 9, 'omitnan');
        bin = discretize(ratio, edges);
        agree = M.match_ref(r, A.P(i).refs{2}(:), Fs, tol) & M.match_ref(r, A.P(i).refs{3}(:), Fs, tol);
        oU = mu & ~mm; oM = mm & ~mu;
        for b = 1:nb
            cntU(b) = cntU(b) + sum(oU & bin == b);
            cntM(b) = cntM(b) + sum(oM & bin == b);
            tot(b)  = tot(b)  + sum(bin == b);
        end
        agU = agU + [sum(oU & agree), sum(oU & ~agree)];
        agM = agM + [sum(oM & agree), sum(oM & ~agree)];
        if sum(agree) >= 20
            sU(end+1) = sum(mu & agree) / sum(agree); %#ok<AGROW>
            sM(end+1) = sum(mm & agree) / sum(agree); %#ok<AGROW>
        end
    end
    fprintf('\n=== MIMIC AF database, %s ===\n', gn{grp});
    fprintf('%-14s %8s %8s %8s\n', 'RR/median RR', 'ref', 'USER-only', 'PB-only');
    for b = 1:nb
        fprintf('%5.1f-%-8.1f %8d %8d %8d\n', edges(b), edges(b+1), tot(b), cntU(b), cntM(b));
    end
    fprintf('Exclusive beats where all 3 refs agree:   USER %d  PB %d\n', agU(1), agM(1));
    fprintf('Exclusive beats where refs disagree:      USER %d  PB %d\n', agU(2), agM(2));
    if numel(sU) >= 5
        fprintf('Sensitivity on agreed beats: %.3f vs %.3f  (p = %.3g, n = %d recordings)\n', ...
            mean(sU), mean(sM), signrank(sU, sM), numel(sU));
    end
end
end

% =========================================================================
% LoA. HEART-RATE LIMITS OF AGREEMENT ON THE HELD-OUT SET
% =========================================================================
function loa_check(ACC, cfg)
% Both detectors against one common reference heart rate, on recordings
% where both return a heart rate, stratified by reference agreement.
M = metrics();
k = find(strcmp({ACC.name}, 'mimic_large_test'), 1);
if isempty(k), fprintf('mimic_large_test not run, skipped.\n'); return; end
A = ACC(k); Fs = A.Fs; P = A.P;
cfg.LAG_WIN_SEC = M.lag_window(A.name);
du = find(strcmpi(A.dets, 'user'), 1); dm = find(strcmpi(A.dets, 'msptd_pb'), 1);
n = numel(P);
eU = nan(n,1); eM = nan(n,1); eE = nan(n,1); md = nan(n,1); vU = false(n,1); vM = false(n,1);
for i = 1:n
    r = P(i).refs; c = cellfun(@numel, r);
    if min(c) < 10, continue; end
    md(i) = (max(c) - min(c)) / max(c);
    hr_common = 60 / mean(M.clean_rr(diff(r{1}) / Fs));
    mu = M.compute(r{1}, P(i).dets{du}, Fs, A.tol, 0, P(i).dur, cfg);
    mm = M.compute(r{1}, P(i).dets{dm}, Fs, A.tol, 0, P(i).dur, cfg);
    me = M.compute(r{1}, r{2},          Fs, A.tol, 0, P(i).dur, cfg);   % PT vs MinEl ECG
    vU(i) = ~isnan(mu.hr_ppg); vM(i) = ~isnan(mm.hr_ppg);
    if vU(i) && vM(i)
        eU(i) = hr_common - mu.hr_ppg;
        eM(i) = hr_common - mm.hr_ppg;
    end
    if ~isnan(me.hr_ppg), eE(i) = hr_common - me.hr_ppg; end
end
fprintf('Valid heart rate: USER only %d | PB only %d | both %d\n', sum(vU & ~vM), sum(vM & ~vU), sum(vU & vM));
loa = @(e) sprintf('bias %+.2f  LoA [%.2f, %.2f]  n=%d', mean(e, 'omitnan'), ...
    mean(e, 'omitnan') - 1.96*std(e, 'omitnan'), mean(e, 'omitnan') + 1.96*std(e, 'omitnan'), sum(~isnan(e)));
for thr = [Inf 0.10 0.05]
    keep = md <= thr & ~isnan(eU) & ~isnan(eM);
    fprintf('\nRecordings with reference disagreement <= %g (both detectors valid):\n', thr);
    fprintf('  MinEl_PPG vs PT : %s\n', loa(eU(keep)));
    fprintf('  MSPTDfast vs PT : %s\n', loa(eM(keep)));
    fprintf('  MinEl ECG vs PT : %s   (reference noise floor)\n', loa(eE(keep)));
end
end
