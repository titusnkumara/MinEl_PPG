/* minel_ppg.c - C port of minel_ppg.m (MinEl_PPG beat detector).
 *
 *   int minel_ppg(const double *ppg, int n, double fs, int *beats);
 *   int minel_ppg_fs(const double *ppg, int n, double fs, double fs_min, int *beats);
 *
 * fs_min: lowest rate kept by the decimation step (every D-th sample, D the
 *         largest integer with fs/D >= fs_min); minel_ppg() uses 30 Hz, as
 *         minel_ppg.m does by default. Values below 30 Hz are rejected (-2).
 * ppg   : raw PPG, n samples at fs Hz
 * beats : output, native 1-based sample indices (MATLAB convention);
 *         must hold at least n/2 + 1 entries
 * returns the number of beats (or -1 on allocation failure)
 *
 * The algorithm and every constant are those of minel_ppg.m; see that file
 * for the description. Internal arrays use 1-based indexing so each line can
 * be checked against the MATLAB source.
 */
#include <stdlib.h>
#include <string.h>
#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846   /* not defined in strict C/POSIX modes (e.g. macOS) */
#endif
#include "numerics.h"
#include "minel_ppg.h"

#define MAXV(a,b) ((a) > (b) ? (a) : (b))
#define MINV(a,b) ((a) < (b) ? (a) : (b))

/* ---------------------------------------------------------------------- */
/* 0.5-5 Hz first-order Butterworth band-pass, cached per sampling rate.
 * Computed along the same path as butter(1, [0.5 5]/(fs/2), 'bandpass'):
 * analog prototype pole -1, band-pass transform (sftrans), bilinear
 * transform with T = 2, polynomial expansion (poly). The cache is thread-local,
 * so the detector can run in several threads at once. Following this path,
 * rather than a closed-form expression, keeps the coefficients identical to
 * the last bit, which matters: a 1e-16 difference can change which of two
 * near-equal samples is a block extremum. */
#if defined(__STDC_VERSION__) && __STDC_VERSION__ >= 201112L
#define TLS _Thread_local
#else
#define TLS __thread
#endif
static void bandpass_coeffs(double fs, double b[3], double a[3])
{
    static TLS double c_fs = -1, cb[3], ca[3];         /* one cache per thread */
    if (fs != c_fs) {
        const double T = 2.0;
        double w1 = 0.5 / (fs / 2), w2 = 5.0 / (fs / 2);
        double Fl = 2 / T * tan(M_PI * w1 / T), Fh = 2 / T * tan(M_PI * w2 / T);
        /* sftrans (band-pass): one pole at -1, no zeros */
        double Sg = pow(1.0 / (Fh - Fl), -1.0);
        double bb = -1.0 * ((Fh - Fl) / 2.0);
        double r  = sqrt(bb * bb - Fh * Fl);
        double Sp1 = bb + r, Sp2 = bb - r;               /* zero at s = 0 */
        /* bilinear, T = 2 */
        double Zg = Sg * ((2 - 0.0 * T) / T) / (((2 - Sp1 * T) / T) * ((2 - Sp2 * T) / T));
        double Zp1 = (2 + Sp1 * T) / (2 - Sp1 * T), Zp2 = (2 + Sp2 * T) / (2 - Sp2 * T);
        /* zeros at z = 1 and z = -1: poly([1 -1]) = [1 0 -1] */
        cb[0] = Zg * 1.0; cb[1] = Zg * 0.0; cb[2] = Zg * -1.0;
        /* poly([Zp1 Zp2]) */
        ca[0] = 1.0; ca[1] = (0.0 - Zp1 * 1.0) - Zp2 * 1.0; ca[2] = 0.0 - Zp2 * (0.0 - Zp1 * 1.0);
        c_fs = fs;
    }
    memcpy(b, cb, 3 * sizeof(double)); memcpy(a, ca, 3 * sizeof(double));
}

/* public: the same band-pass coefficients, for use by other pipelines */
void minel_bandpass(double fs, double b[3], double a[3]) { bandpass_coeffs(fs, b, a); }

/* largest integer stride D that keeps fs/D >= fs_min (1 if none) */
static int decimation_stride(double fs, double fs_min)
{
    double k = floor(fs / fs_min);
    return k >= 1 ? (int)k : 1;
}

/* ---------------------------------------------------------------------- */
/* position (1-based) of value v in sorted array s of length n, first match; 0 if absent */
static int find_pos(const int *s, int n, int v)
{
    int lo = 0, hi = n - 1;
    while (lo < hi) { int mid = (lo + hi) / 2; if (s[mid] < v) lo = mid + 1; else hi = mid; }
    return (n > 0 && s[lo] == v) ? lo + 1 : 0;
}

/* ---------------------------------------------------------------------- */
/* rhythm-aware detector on the sparse stream; returns 1-based positions */
static int detect_beats(const double *s1, int N, double fs, int *out)
{
    if (N < 15) return 0;
    const double AMP_CAP = 1.5, AMP_W = 0.3, RR_G1 = 0.6, RR_G2 = 1.6;
    const int GATE_ESCAPE = 3; const double HR_INIT = 120, HR_UP = 130, HR_DN = 110;
    int count = 0;

    double *s = malloc(sizeof(double) * (N + 2)), *sa = malloc(sizeof(double) * (N + 3));
    for (int i = 1; i <= N; i++) { s[i] = s1[i-1]; sa[i] = fabs(s[i]); }
    sa[0] = 0; sa[N+1] = 0; sa[N+2] = 0;

    double msum = 0; for (int i = 1; i <= N; i++) msum += sa[i];
    double off = 0.0010 * (msum / N), alpha = 0.202;
    char *m_slow = calloc(N + 2, 1), *m_fast = calloc(N + 2, 1), *fast = calloc(N + 2, 1), *mask = calloc(N + 2, 1);
    for (int i = 1; i <= N; i++) {
        /* conv(x, ones(1,k)/k, 'same'), summed in filter() order (oldest tap innermost) */
        double u1 = (i + 1 <= N) ? sa[i+1] : 0.0;
        double c2 = 0.5 * sa[i] + 0.5 * u1;
        double um2 = (i - 2 >= 1) ? sa[i-2] : 0.0, um1 = (i - 1 >= 1) ? sa[i-1] : 0.0;
        double up1 = (i + 1 <= N) ? sa[i+1] : 0.0, up2 = (i + 2 <= N) ? sa[i+2] : 0.0;
        const double f5 = 1.0 / 5.0, f3 = 1.0 / 3.0;
        double c5 = (((f5 * um2 + f5 * um1) + f5 * sa[i]) + f5 * up1) + f5 * up2;
        double c3 = (f3 * um1 + f3 * sa[i]) + f3 * up1;
        m_slow[i] = (alpha * sa[i] + (1 - alpha) * c2) > c5 + off;
        m_fast[i] = sa[i] > c3 + off;
    }

    /* heart-rate regime with hysteresis */
    int *st = malloc(sizeof(int) * (N + 1)), nst = 0;
    for (int i = 1; i <= N; i++) if (m_slow[i] && (i == 1 || !m_slow[i-1])) st[++nst] = i;
    if (nst >= 4) {
        int nr = nst - 1;
        double *dst = malloc(sizeof(double) * nr), *mm = malloc(sizeof(double) * nr);
        char *reg = malloc(nr + 1);
        for (int k = 0; k < nr; k++) dst[k] = st[k+2] - st[k+1];
        movmedian_odd(dst, nr, 9, mm);
        int state = (60 * fs / mm[0]) > HR_INIT;
        for (int k = 1; k <= nr; k++) {
            double hr = 60 * fs / mm[k-1];
            if (state && hr < HR_DN) state = 0; else if (!state && hr > HR_UP) state = 1;
            reg[k] = (char)state;
        }
        int cs = 0, nxt = 2;                  /* cs = cumsum(mk), mk(st(2:end)) = true */
        for (int i = 1; i <= N; i++) {
            if (nxt <= nst && st[nxt] == i) { cs++; nxt++; }
            fast[i] = reg[MINV(cs + 1, nr)];
        }
        free(dst); free(mm); free(reg);
    }
    for (int i = 1; i <= N; i++) mask[i] = fast[i] ? m_fast[i] : m_slow[i];

    /* candidates: maximum of each region */
    int *cl = malloc(sizeof(int) * (N + 1)); double *ca = malloc(sizeof(double) * (N + 1)); int nc = 0;
    for (int i = 1; i <= N; i++) {
        if (mask[i] && (i == 1 || !mask[i-1])) {
            int e = i; while (e + 1 <= N && mask[e+1]) e++;
            int im = i; for (int k = i + 1; k <= e; k++) if (s[k] > s[im]) im = k;
            nc++; cl[nc] = im; ca[nc] = s[im];
        }
    }
    if (nc == 0) goto done;

    int win_len = MINV(N, (int)round(58 * fs));
    double base = 0.001 * prctile_dbl(s + 1 + (N - win_len), win_len, 71);
    int nv = 0;
    for (int k = 1; k <= nc; k++) if (ca[k] > base) { nv++; cl[nv] = cl[k]; ca[nv] = ca[k]; }
    nc = nv;
    if (nc == 0) goto done;

    /* RR-adaptive acceptance with the dicrotic contest */
    const int min_dist = 1; const double rr_tol = 0.123, low_f = 0.366, high_f = 1.220;
    const double MIN_RR_FLOOR = 0.35 * fs;
    int cap_q = nc + 4;
    int *qi = malloc(sizeof(int) * (cap_q + 1)); double *qa = malloc(sizeof(double) * (cap_q + 1));
    int acc = 0;
    double rr_hist[11] = {0}, amp_hist[11] = {0}; int rr_cnt = 0, rr_p = 0, amp_cnt = 0, amp_p = 0;
    int last_rr_added = 0, n_gate_miss = 0; double dev = 0, med_rr = 0;
    int rr_dirty = 1, cv_dirty = 1, amp_dirty = 1; double rr_med = 0, cv_rr = 0, amp_med = 0;

    for (int i = 1; i <= nc; i++) {
        int cur_t = cl[i]; double cur_a = ca[i];
        int cf = fast[cur_t];
        int have_rr = acc >= 2 && rr_cnt >= 1;
        double relax = 0;
        if (have_rr && rr_dirty) { rr_med = median_dbl(rr_hist + 1, rr_cnt); rr_dirty = 0; }
        if (have_rr && rr_cnt >= 5 && !cf) {
            if (cv_dirty) {
                double dvv[10];
                for (int k = 0; k < rr_cnt; k++) dvv[k] = fabs(rr_hist[k+1] - rr_med);
                double rr_mad = 1.4826 * median_dbl(dvv, rr_cnt);
                cv_rr = rr_mad / MAXV(rr_med, 2.220446049250313e-16);
                cv_dirty = 0;
            }
            relax = MINV(1.0, MAXV(0.0, (cv_rr - 0.10) / 0.15));
        }
        double notch_thr = 0.60 - 0.25 * relax;
        double rr_dev = NAN;
        if (have_rr) { med_rr = rr_med; rr_dev = fabs(cur_t - (qi[acc] + med_rr)) / med_rr; }
        double mult = low_f; if (!isnan(rr_dev) && rr_dev > rr_tol) mult = high_f;
        if (acc >= 1 && (cur_t - qi[acc]) < min_dist) continue;
        if (cur_a <= mult * base) continue;

        if (have_rr) {
            int dt = cur_t - qi[acc];
            double lo = cf ? MAXV(0.35 * med_rr, MIN_RR_FLOOR) : MAXV(0.30 * med_rr, MIN_RR_FLOOR);
            double hi = 0.70 * med_rr;
            if (dt >= lo && dt <= hi) {
                if (amp_dirty) { amp_med = median_dbl(amp_hist + 1, amp_cnt); amp_dirty = 0; }
                double ref_a = amp_med;
                int prev_t = qi[acc-1];
                double cv = 1.25 * dev / med_rr;
                double tw = MINV(1.0, MAXV(0.0, (med_rr - 4) / 4)) * MINV(1.0, MAXV(0.0, (0.25 - cv) / 0.15));
                double w = AMP_W * tw + (1 - tw), cap = AMP_CAP * tw + 10 * (1 - tw);
                double sg = MAXV(1.0, 0.2 * med_rr);
                double dl = (qi[acc] - prev_t - med_rr) / sg, dc = (cur_t - prev_t - med_rr) / sg;
                double s_last = pow(exp(-0.5 * dl * dl), tw) * pow(MINV(MAXV(qa[acc] / ref_a, 0.05), cap), w);
                double s_cur  = pow(exp(-0.5 * dc * dc), tw) * pow(MINV(MAXV(cur_a / ref_a, 0.05), cap), w);
                if (s_cur > s_last) {
                    qi[acc] = cur_t; qa[acc] = cur_a;
                    if (last_rr_added) { rr_hist[rr_p] = cur_t - prev_t; rr_dirty = cv_dirty = 1; }
                    amp_hist[amp_p] = cur_a; amp_dirty = 1;
                    continue;
                }
                if (cur_a / ref_a < notch_thr) continue;
            }
        }
        /* accept */
        acc++; qi[acc] = cur_t; qa[acc] = cur_a;
        amp_p = amp_p % 10 + 1; amp_hist[amp_p] = cur_a; amp_cnt = MINV(amp_cnt + 1, 10);
        amp_dirty = 1;
        last_rr_added = 0;
        if (acc >= 2) {
            double new_rr = qi[acc] - qi[acc-1];
            int ok = 1;
            if (rr_cnt >= 3) {
                ok = new_rr >= RR_G1 * med_rr && new_rr <= RR_G2 * med_rr;
                if (ok) n_gate_miss = 0;
                else { n_gate_miss++; if (n_gate_miss >= GATE_ESCAPE) { rr_cnt = 0; rr_p = 0; n_gate_miss = 0; ok = 1; rr_dirty = cv_dirty = 1; } }
            }
            if (ok) {
                if (rr_cnt == 0) dev = 0.15 * new_rr; else dev = 0.8 * dev + 0.2 * fabs(new_rr - med_rr);
                rr_p = rr_p % 10 + 1; rr_hist[rr_p] = new_rr; rr_cnt = MINV(rr_cnt + 1, 10);
                last_rr_added = 1; rr_dirty = cv_dirty = 1;
            }
        }
    }

    /* searchback */
    if (acc >= 2) {
        double *rr_all = malloc(sizeof(double) * (cap_q + 1));
        int nra = acc - 1;
        for (int k = 1; k <= nra; k++) rr_all[k] = qi[k+1] - qi[k];
        int s0 = MAXV(1, nra - 9);
        double mrr = median_dbl(rr_all + s0, nra - s0 + 1);
        double ref_a = median_dbl(qa + 1, acc);
        double m_all = median_dbl(rr_all + 1, nra), *dv = malloc(sizeof(double) * (nra + 1));
        for (int k = 1; k <= nra; k++) dv[k] = fabs(rr_all[k] - m_all);
        double cv_g = 1.4826 * median_dbl(dv + 1, nra) / MAXV(m_all, 2.220446049250313e-16);
        free(dv);
        int idx = 2;
        while (idx <= acc) {
            int a = qi[idx-1], b = qi[idx];
            double mlt = fast[a] ? 1.8 : 1.5;
            if (b - a > mlt * mrr) {
                double tw = MINV(1.0, MAXV(0.0, (mrr - 4) / 4)) * MINV(1.0, MAXV(0.0, (0.25 - cv_g) / 0.15));
                double w = AMP_W * tw + (1 - tw), cap = AMP_CAP * tw + 10 * (1 - tw), sg = MAXV(1.0, 0.2 * mrr);
                int best = 0; double bsc = -INFINITY;
                for (int k = 1; k <= nc; k++) {
                    if (cl[k] >= a + min_dist && cl[k] <= b - min_dist) {
                        double z = (cl[k] - a - mrr) / sg;
                        double sc = pow(exp(-0.5 * z * z), tw) * pow(MINV(MAXV(ca[k] / ref_a, 0.05), cap), w);
                        if (best == 0 || sc > bsc) { best = k; bsc = sc; }
                    }
                }
                if (best && ca[best] > 0.5 * base) {
                    if (acc + 1 > cap_q) {
                        cap_q *= 2;
                        qi = realloc(qi, sizeof(int) * (cap_q + 1)); qa = realloc(qa, sizeof(double) * (cap_q + 1));
                        rr_all = realloc(rr_all, sizeof(double) * (cap_q + 1));
                    }
                    memmove(qi + idx + 1, qi + idx, sizeof(int) * (acc - idx + 1));
                    memmove(qa + idx + 1, qa + idx, sizeof(double) * (acc - idx + 1));
                    qi[idx] = cl[best]; qa[idx] = ca[best]; acc++;
                    nra = acc - 1;
                    for (int k = 1; k <= nra; k++) rr_all[k] = qi[k+1] - qi[k];
                    s0 = MAXV(1, nra - 9);
                    mrr = median_dbl(rr_all + s0, nra - s0 + 1);
                    continue;
                }
            }
            idx++;
        }
        free(rr_all);
    }

    /* final close-pair resolution */
    if (acc >= 3) {
        int nra = acc - 1;
        double *rr_all = malloc(sizeof(double) * (nra + 1)), *dv = malloc(sizeof(double) * (nra + 1));
        for (int k = 1; k <= nra; k++) rr_all[k] = qi[k+1] - qi[k];
        double mrr = median_dbl(rr_all + 1, nra);
        int rf = MAXV((int)round(0.25 * fs), (int)round(0.35 * mrr));
        int rs = MAXV((int)round(0.30 * fs), (int)round(0.40 * mrr));
        double ref_a = median_dbl(qa + 1, acc);
        for (int k = 1; k <= nra; k++) dv[k] = fabs(rr_all[k] - mrr);
        double cv_g = 1.4826 * median_dbl(dv + 1, nra) / MAXV(mrr, 2.220446049250313e-16);
        double tw = MINV(1.0, MAXV(0.0, (mrr - 4) / 4)) * MINV(1.0, MAXV(0.0, (0.25 - cv_g) / 0.15));
        double w = AMP_W * tw + (1 - tw), cap = AMP_CAP * tw + 10 * (1 - tw), sg = MAXV(1.0, 0.2 * mrr);
        int nq = acc, i = 2;
        while (i <= nq) {
            int refrac = fast[qi[i]] ? rf : rs;
            if (qi[i] - qi[i-1] < refrac) {
                int drop_second;
                if (i > 2) {
                    int anc = qi[i-2];
                    double d1 = (qi[i-1] - anc - mrr) / sg, d2 = (qi[i] - anc - mrr) / sg;
                    double sc1 = pow(exp(-0.5 * d1 * d1), tw) * pow(MINV(MAXV(qa[i-1] / ref_a, 0.05), cap), w);
                    double sc2 = pow(exp(-0.5 * d2 * d2), tw) * pow(MINV(MAXV(qa[i] / ref_a, 0.05), cap), w);
                    drop_second = sc2 < sc1;
                } else drop_second = qa[i] < qa[i-1];
                int del = drop_second ? i : i - 1;
                memmove(qi + del, qi + del + 1, sizeof(int) * (nq - del));
                memmove(qa + del, qa + del + 1, sizeof(double) * (nq - del));
                nq--;
                if (!drop_second && i > 2) i--;
            } else i++;
        }
        acc = nq;
        free(rr_all); free(dv);
    }
    for (int k = 1; k <= acc; k++) out[k-1] = qi[k];
    count = acc;
    free(qi); free(qa);
done:
    free(s); free(sa); free(m_slow); free(m_fast); free(fast); free(mask); free(st); free(cl); free(ca);
    return count;
}

/* ---------------------------------------------------------------------- */
/* adaptive refractory blanking (cand sorted, 1-based sample indices) */
static int blank_refractory(int *cand, int n, const int *sidx, const double *sval, int ns, double fs)
{
    const double SWAP_RATIO = 1.5;
    int *c = malloc(sizeof(int) * (n + 1)); double *amp = malloc(sizeof(double) * (n + 1));
    int m = 0;
    for (int i = 0; i < n; i++) {
        int p = find_pos(sidx, ns, cand[i]);
        if (p) { m++; c[m] = cand[i]; amp[m] = sval[p-1]; }
    }
    if (m == 0) { free(c); free(amp); return 0; }
    int MIN_RR = MAXV(2, (int)round(0.27 * fs)), MAX_REFR = (int)round(1.00 * fs);
    int MIN_RR_FLOOR = (int)round(0.35 * fs);
    const double ADAPT = 0.45, MUST = 1.10, PLO = 0.40, PHI = 0.65, AMPR = 0.40; const int AMP_MINB = 5, RW = 8;
    double rr_hist[8], amp_hist[8]; int rr_count = 0, amp_count = 1;
    amp_hist[0] = amp[1];
    char *accd = calloc(m + 1, 1); accd[1] = 1;
    int last_i = 1; int prev_peak = -1; int have_prev = 0;
    int rr_dirty = 1, amp_dirty = 1; double rr_med = 0, amp_med = 0;
    for (int i = 2; i <= m; i++) {
        int dt = c[i] - c[last_i];
        if (rr_count > 0 && rr_dirty) { rr_med = median_dbl(rr_hist, rr_count); rr_dirty = 0; }
        double med_rr = rr_count == 0 ? 2.0 * MIN_RR : rr_med;
        int rej = dt < MIN_RR;
        if (!rej && dt <= MUST * med_rr) {
            int refr = MINV(MAX_REFR, MAXV(MIN_RR, (int)round(ADAPT * med_rr)));
            if (dt < refr) rej = 1;
            else if (amp_count >= AMP_MINB) {
                if (amp_dirty) { amp_med = median_dbl(amp_hist, amp_count); amp_dirty = 0; }
                double med_amp = amp_med;
                double lo = MAXV(PLO * med_rr, (double)MIN_RR_FLOOR), hi = PHI * med_rr;
                if (med_amp > 0 && dt >= lo && dt <= hi && amp[i] < AMPR * med_amp) rej = 1;
            }
        }
        if (rej && amp[i] > SWAP_RATIO * amp[last_i]) {
            if (!have_prev || (c[i] - prev_peak) >= MIN_RR) {
                accd[last_i] = 0; accd[i] = 1;
                if (have_prev && rr_count > 0) { rr_hist[rr_count-1] = c[i] - prev_peak; rr_dirty = 1; }
                if (amp_count > 0) { amp_hist[amp_count-1] = amp[i]; amp_dirty = 1; }
                last_i = i;
                continue;
            }
        }
        if (rej) continue;
        accd[i] = 1; prev_peak = c[last_i]; have_prev = 1;
        if (rr_count < RW) rr_hist[rr_count++] = dt;
        else { memmove(rr_hist, rr_hist + 1, sizeof(double) * (RW - 1)); rr_hist[RW-1] = dt; }
        if (amp_count < RW) amp_hist[amp_count++] = amp[i];
        else { memmove(amp_hist, amp_hist + 1, sizeof(double) * (RW - 1)); amp_hist[RW-1] = amp[i]; }
        rr_dirty = amp_dirty = 1;
        last_i = i;
    }
    int k = 0;
    for (int i = 1; i <= m; i++) if (accd[i]) cand[k++] = c[i];
    free(c); free(amp); free(accd);
    return k;
}

/* movmad(x, 30) with a sliding sorted window: each step inserts one value
 * and removes one, then reads the median and the median absolute deviation
 * by merging outwards from the median. Same multiset and same median
 * formula, so the result equals movmad_even() exactly. */
static void movmad30_sliding(const double *x, int n, double *out)
{
    const int kb = 15, kf = 14;
    double s[32]; int m = 0;
    for (int i = 0; i <= kf && i < n; i++) {                 /* initial window [0, kf] */
        double v = x[i]; int j = m - 1;
        while (j >= 0 && s[j] > v) { s[j+1] = s[j]; j--; }
        s[j+1] = v; m++;
    }
    for (int i = 0; i < n; i++) {
        if (i > 0) {
            int add = i + kf, rem = i - kb - 1;
            if (rem >= 0) {                                   /* remove x[rem] */
                double v = x[rem]; int j = 0;
                while (s[j] != v) j++;
                for (; j < m - 1; j++) s[j] = s[j+1];
                m--;
            }
            if (add < n) {                                    /* insert x[add] */
                double v = x[add]; int j = m - 1;
                while (j >= 0 && s[j] > v) { s[j+1] = s[j]; j--; }
                s[j+1] = v; m++;
            }
        }
        double med = (m % 2) ? s[m/2] : (s[m/2 - 1] + s[m/2]) / 2.0;
        int p = 0; while (p < m && s[p] < med) p++;           /* first element >= median */
        int l = p - 1, r = p, need = m / 2;
        double prev = 0, cur = 0;
        for (int k = 0; k <= need; k++) {                      /* k-th smallest deviation */
            double dl = (l >= 0) ? fabs(s[l] - med) : INFINITY;
            double dr = (r < m)  ? fabs(s[r] - med) : INFINITY;
            prev = cur;
            if (dl <= dr) { cur = dl; l--; } else { cur = dr; r++; }
        }
        out[i] = (m % 2) ? cur : (prev + cur) / 2.0;
    }
}

/* slope and amplitude quality gate */
static int quality_gate(int *cand, int n, const int *sidx, const double *sval, int ns, double fs)
{
    if (n == 0) return 0;
    double *slopes = malloc(sizeof(double) * n), *mad = malloc(sizeof(double) * n);
    int *pos = malloc(sizeof(int) * n);
    for (int k = 0; k < n; k++) {
        int p = find_pos(sidx, ns, cand[k]); pos[k] = p;
        if (p > 1) {
            double amp = sval[p-1] - sval[p-2];
            double dt = (sidx[p-1] - sidx[p-2]) / fs;
            slopes[k] = MAXV(0.0, amp / dt);
        } else {
            slopes[k] = fabs(sval[p-1] - sval[p]) / ((sidx[p] - sidx[p-1]) / fs);
        }
    }
    double sqi;
    int m = 0;
    if (n >= 30) {
        movmad30_sliding(slopes, n, mad);
        for (int k = 0; k < n; k++) {
            double md = mad[k] < 1e-6 ? 1e-6 : mad[k];
            sqi = MINV(1.0, slopes[k] / (4 * md));
            if (sqi >= 0.025) { cand[m] = cand[k]; pos[m] = pos[k]; m++; }
        }
    } else {
        double med = median_dbl(slopes, n), *dv = malloc(sizeof(double) * n);
        for (int k = 0; k < n; k++) dv[k] = fabs(slopes[k] - med);
        double md = median_dbl(dv, n); if (md < 1e-6) md = 1e-6;
        free(dv);
        for (int k = 0; k < n; k++) {
            sqi = MINV(1.0, slopes[k] / (4 * md));
            if (sqi >= 0.025) { cand[m] = cand[k]; pos[m] = pos[k]; m++; }
        }
    }
    if (m > 0) {
        double *h = malloc(sizeof(double) * m);
        for (int k = 0; k < m; k++) h[k] = sval[pos[k]-1];
        double thr = 0.025 * prctile_dbl(h, m, 70);
        int q = 0;
        for (int k = 0; k < m; k++) if (h[k] >= thr) cand[q++] = cand[k];
        m = q; free(h);
    }
    free(slopes); free(mad); free(pos);
    return m;
}

/* ---------------------------------------------------------------------- */
int minel_ppg(const double *ppg, int n, double fs, int *beats)
{
    return minel_ppg_fs(ppg, n, fs, MINEL_FS_MIN_DEFAULT, beats);
}

int minel_ppg_fs(const double *ppg, int n, double fs, double fs_min, int *beats)
{
    if (!(fs_min >= MINEL_FS_MIN_DEFAULT)) return -2;
    int D = decimation_stride(fs, fs_min);
    double fsu = fs / D;
    int nd = (n + D - 1) / D;
    double *dx = malloc(sizeof(double) * (nd + 1)), *d = malloc(sizeof(double) * (nd + 1));
    if (!dx || !d) return -1;
    double prev = isfinite(ppg[0]) ? ppg[0] : 0.0;           /* decimate + difference, one pass */
    dx[0] = 0;
    for (int i = 1; i < nd; i++) {
        double v = ppg[(long)i * D]; if (!isfinite(v)) v = 0.0;
        dx[i] = v - prev; prev = v;
    }
    double b[3], a[3];
    bandpass_coeffs(fsu, b, a);
    int nb_out = 0;

    /* 0.2 s block boundaries, same arithmetic as sparse_extrema():
       integer rates in integer arithmetic, non-integer rates (e.g. 31.25 Hz
       after decimation) as round(j * 2*fs/T) */
    const int T = 10;
    int integer_rate = (fsu == round(fsu)), FsI = (int)round(fsu);
    int nblk = integer_rate ? (int)floor((double)nd * T / (2.0 * FsI))
                            : (int)floor((double)nd * T / (2.0 * fsu));
    int *sidx = NULL; double *sval = NULL; int ns = 0;
    if (nblk < 20) goto out;
    long *e = malloc(sizeof(long) * (size_t)(nblk + 1));
    int *imn = malloc(sizeof(int) * (size_t)nblk), *imx = malloc(sizeof(int) * (size_t)nblk);
    if (!e || !imn || !imx) { free(e); free(imn); free(imx); nb_out = -1; goto out; }
    const double cblk = 2.0 * fsu / T;                       /* block length in samples */
    for (int j = 0; j <= nblk; j++) {
        if (integer_rate)
            e[j] = ((2 * FsI) % T == 0) ? (long)j * (2 * FsI / T) : (long)floor((4.0 * j * FsI + T) / (2.0 * T));
        else
            e[j] = (long)round(j * cblk);
    }
    /* band-pass filtfilt with the block extrema found during the reverse pass */
    if (filtfilt3_extrema(b, a, dx, nd, d, e, nblk, imn, imx) != 0) { free(e); free(imn); free(imx); goto out; }
    ns = 2 * nblk;
    sidx = malloc(sizeof(int) * (size_t)ns); sval = malloc(sizeof(double) * (size_t)ns);
    for (int j = 0; j < nblk; j++) {                          /* blocks are ordered: no sort needed */
        int lo = MINV(imn[j], imx[j]), hi = MAXV(imn[j], imx[j]);
        sidx[2*j] = lo; sidx[2*j+1] = hi;
        sval[2*j] = d[lo-1]; sval[2*j+1] = d[hi-1];
    }
    free(e); free(imn); free(imx);

    int *pos = malloc(sizeof(int) * (ns + 1));
    int np = detect_beats(sval, ns, 10.0, pos);
    int *cand = malloc(sizeof(int) * (np + 1));
    for (int k = 0; k < np; k++) cand[k] = sidx[pos[k] - 1];
    sort_int(cand, np);
    np = blank_refractory(cand, np, sidx, sval, ns, fsu);
    np = quality_gate(cand, np, sidx, sval, ns, fsu);
    for (int k = 0; k < np; k++) {                       /* parabolic refinement */
        double pkf = cand[k];
        int kk = cand[k];
        if (kk > 1 && kk < nd) {
            double ym = d[kk-2], y0 = d[kk-1], yp = d[kk];
            double den = ym - 2 * y0 + yp, delta = 0;
            if (den < 0) delta = 0.5 * (ym - yp) / den;
            pkf = kk + MAXV(-0.5, MINV(0.5, delta));
        }
        beats[k] = (int)round(D * (pkf - 1) + 1);
    }
    nb_out = unique_int(beats, np);
    free(pos); free(cand);
out:
    free(sidx); free(sval); free(dx); free(d);
    return nb_out;
}
