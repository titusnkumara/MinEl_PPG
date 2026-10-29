/* msptd_pb.c - C port of MSPTDfast v2 as the ppg-beats benchmark runs it
 *
 * LICENSE: GPL-3.0 (not MIT like the rest of this repository).
 * This file translates code from the ppg-beats toolbox by Peter H. Charlton:
 *   msptdfastv2_beat_detector.m (MIT, Copyright (c) 2022, 2024 Peter H. Charlton)
 *   detect_ppg_beats.m, tidy_beats.m (GPL-3.0, Copyright (c) 2022 Peter H. Charlton)
 * As a translation of GPL-3.0 code, this file is distributed under the GNU
 * General Public License v3.0 or later: <https://www.gnu.org/licenses/>.
 * It is provided WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * (Charlton et al 2025; ppg-beats toolbox, Charlton 2024):
 *
 *   1. pre-filter (preprocess_ppg_signal.m): Kaiser-window FIR low-pass
 *      (stop 6.93 Hz, pass 10 Hz), then baseline removal by subtracting a
 *      low-pass (stop 0.6 Hz, pass 1.02 Hz) of the Tukey-windowed signal,
 *      both zero-phase; designed as MATLAB kaiserord + fir1(..., 'scale')
 *   2. msptdfastv2_beat_detector.m: 6 s windows (20% overlap), decimation
 *      to about 20 Hz, linear detrend, multi-scale local maxima/minima over
 *      scales covering heart rates >= 30 bpm, refinement at the native rate
 *   3. detect_ppg_beats.m: tidy_peaks_and_onsets and calc_mid_amp_points
 *      (the mid-amplitude fiducial the benchmark scores)
 *
 * Each function names the toolbox function it reproduces. Indices are
 * 1-based sample indices (MATLAB convention) throughout.
 */
#include <stdlib.h>
#include <string.h>
#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846   /* not defined in strict C/POSIX modes (e.g. macOS) */
#endif
#include "numerics.h"
#include "msptd_pb.h"

#define MAXV(a,b) ((a) > (b) ? (a) : (b))
#define MINV(a,b) ((a) < (b) ? (a) : (b))

typedef struct { int *v; int n, cap; } ivec;
static void iv_push(ivec *a, int x)
{
    if (a->n == a->cap) { a->cap = a->cap ? 2 * a->cap : 64; a->v = realloc(a->v, sizeof(int) * (size_t)a->cap); }
    a->v[a->n++] = x;
}

/* ====================================================================== */
/* 1. PRE-FILTER                                                          */
/* ====================================================================== */
static double bessel_i0(double x)
{
    double sum = 1.0, term = 1.0, q = x * x / 4.0;
    for (int k = 1; k < 200; k++) {
        term *= q / ((double)k * k);
        sum += term;
        if (term < 1e-17 * sum) break;
    }
    return sum;
}

/* MATLAB: [N,Wn,beta,type] = kaiserord([Fstop Fpass]/(fs/2), [1 0], [0.01 0.05]);
 *         b = fir1(N, Wn, type, kaiser(N+1, beta), 'scale');            */
static double *kaiser_lpf(double Fstop, double Fpass, double fs, int *ntaps)
{
    double A = -20.0 * log10(0.01);                      /* min deviation */
    double beta = (A > 50) ? 0.1102 * (A - 8.7)
                : (A >= 21) ? 0.5842 * pow(A - 21, 0.4) + 0.07886 * (A - 21) : 0.0;
    double df = (Fpass - Fstop) / fs;                     /* transition width, cycles/sample */
    int N = (int)ceil((A - 7.95) / (2 * M_PI * 2.285 * df));
    double Wn = (Fstop + Fpass) / fs;                     /* cutoff, normalised to Nyquist */
    double *b = malloc(sizeof(double) * (size_t)(N + 1)), s = 0, i0b = bessel_i0(beta);
    for (int k = 0; k <= N; k++) {
        double m = k - N / 2.0;
        double h = (m == 0) ? Wn : sin(M_PI * Wn * m) / (M_PI * m);      /* Wn*sinc(Wn*m) */
        double r = (N == 0) ? 0 : 2.0 * k / N - 1.0;
        double w = bessel_i0(beta * sqrt(MAXV(0.0, 1.0 - r * r))) / i0b;
        b[k] = h * w; s += b[k];
    }
    for (int k = 0; k <= N; k++) b[k] /= s;               /* 'scale': unit gain at DC */
    *ntaps = N + 1;
    return b;
}

/* Filter-design cache, one per thread (thread-local), holding up to
 * NSLOT sampling rates so that a thread alternating between databases does
 * not redesign the filters. Designs are identical to an uncached design;
 * only the reuse changes, so outputs are unaffected. */
#if defined(__STDC_VERSION__) && __STDC_VERSION__ >= 201112L
#define TLS _Thread_local
#else
#define TLS __thread
#endif
#define NSLOT 8
typedef struct { double fs; double *bh, *bl; int nh, nl; } fslot;
static TLS fslot c_slot[NSLOT];
static TLS int c_nslot = 0, c_next = 0;
static TLS double c_fs = -1, *c_bh = NULL, *c_bl = NULL;
static TLS int c_nh, c_nl;
static void design(double fs)
{
    if (fs == c_fs) return;
    int k;
    for (k = 0; k < c_nslot; k++) if (c_slot[k].fs == fs) break;
    if (k == c_nslot) {                              /* not cached: design it */
        if (c_nslot < NSLOT) k = c_nslot++;
        else { k = c_next; c_next = (c_next + 1) % NSLOT; free(c_slot[k].bh); free(c_slot[k].bl); }
        c_slot[k].bh = kaiser_lpf(6.93, 10.0, fs, &c_slot[k].nh);
        c_slot[k].bl = kaiser_lpf(0.6, 1.02, fs, &c_slot[k].nl);
        c_slot[k].fs = fs;
    }
    c_bh = c_slot[k].bh; c_nh = c_slot[k].nh; c_bl = c_slot[k].bl; c_nl = c_slot[k].nl; c_fs = fs;
}
void msptd_release_cache(void)                       /* frees the calling thread's cache */
{
    for (int k = 0; k < c_nslot; k++) { free(c_slot[k].bh); free(c_slot[k].bl); }
    c_nslot = 0; c_next = 0; c_fs = -1; c_bh = c_bl = NULL;
}
void msptd_filter_taps(double fs, const double **bh, int *nh, const double **bl, int *nl)
{
    design(fs); *bh = c_bh; *nh = c_nh; *bl = c_bl; *nl = c_nl;
}

/* preprocess_ppg_signal.m: elim_vhfs3 + elim_vlfs */
static void charlton_filter(const double *x, int n, double fs, double *y)
{
    design(fs);
    double *t = malloc(sizeof(double) * (size_t)n), *u = malloc(sizeof(double) * (size_t)n);
    memcpy(y, x, sizeof(double) * (size_t)n);
    if (n > 3 * (c_nh - 1)) { filtfilt(c_bh, NULL, c_nh, y, n, t); memcpy(y, t, sizeof(double) * (size_t)n); }
    if (n > 3 * (c_nl - 1)) {
        double r = 0.5 * 2 / ((n - 1) / fs);              /* tukeywin(n, r): 0.5 s tapers */
        if (r > 1) r = 1;
        int tl = (int)floor(r * (n - 1) / 2) + 1;          /* taper length */
        for (int i = 0; i < n; i++) {
            double w = 1.0;
            if (r > 0) {
                int j = (i < n - 1 - i) ? i : n - 1 - i;    /* mirror */
                if (j < tl) { double tt = j * (1.0 / (n - 1));    /* linspace(0,1,n) */ w = (1 + cos(M_PI * (2 * tt / r - 1))) / 2; }
            }
            u[i] = y[i] * w;
        }
        filtfilt(c_bl, NULL, c_nl, u, n, t);
        for (int i = 0; i < n; i++) y[i] = y[i] - t[i];
    }
    free(t); free(u);
}

/* ====================================================================== */
/* 2. MSPTDfast v2                                                        */
/* ====================================================================== */
/* detect_peaks_and_onsets_using_msptd (lms_calc_method 1, reduced scales) */
static void msptd_scalogram(const double *x0, int N, double fs, ivec *p, ivec *t)
{
    p->n = 0; t->n = 0;
    int L = (int)ceil(N / 2.0) - 1;
    if (L < 1) return;
    double durn = N / fs;
    int first = 0, max_scale = 0;
    for (int k = 1; k <= L; k++)
        if ((double)L / k / durn >= 30.0 / 60.0) { if (!first) first = k; max_scale = k; }
    if (!max_scale) return;

    /* detrend: remove the least-squares straight line */
    double *x = malloc(sizeof(double) * (size_t)(N + 1));
    double tb = (N + 1) / 2.0, xb = 0, sxt = 0, stt = 0;
    for (int i = 1; i <= N; i++) xb += x0[i-1];
    xb /= N;
    for (int i = 1; i <= N; i++) { sxt += (i - tb) * (x0[i-1] - xb); stt += (i - tb) * (i - tb); }
    double slope = sxt / stt;
    for (int i = 1; i <= N; i++) x[i] = x0[i-1] - (xb + slope * (i - tb));

    /* gamma: number of local maxima (minima) at each scale */
    int lam_max = 1, lam_min = 1, gmax_best = -1, gmin_best = -1;
    for (int k = 1; k <= max_scale; k++) {
        int gmx = 0, gmn = 0;
        for (int q = k + 1; q <= N - k; q++) {
            if (x[q] > x[q-k] && x[q] > x[q+k]) gmx++;
            if (x[q] < x[q-k] && x[q] < x[q+k]) gmn++;
        }
        if (gmx > gmax_best) { gmax_best = gmx; lam_max = k; }
        if (gmn > gmin_best) { gmin_best = gmn; lam_min = k; }
    }
    /* points that are local maxima (minima) at every scale first..lambda */
    for (int q = 1; q <= N; q++) {
        int okp = lam_max >= first, okt = lam_min >= first;
        for (int k = first; k <= lam_max && okp; k++)
            okp = (q - k >= 1 && q + k <= N) && x[q] > x[q-k] && x[q] > x[q+k];
        for (int k = first; k <= lam_min && okt; k++)
            okt = (q - k >= 1 && q + k <= N) && x[q] < x[q-k] && x[q] < x[q+k];
        if (okp) iv_push(p, q);
        if (okt) iv_push(t, q);
    }
    free(x);
}

static TLS int n_bound_guards = 0;          /* refinement windows that MATLAB would reject */

/* msptdpcref_beat_detector with the msptdfastv2 options */
static void msptdfastv2(const double *sig, int n, double fs, ivec *peaks, ivec *onsets)
{
    peaks->n = 0; onsets->n = 0;
    double no = 6.0 * fs;                                  /* samples per window */
    ivec ws = {0}, we = {0};
    if (n <= no) { iv_push(&ws, 1); iv_push(&we, n); }
    else {
        int off = (int)round(no * (1 - 0.2));
        for (double s = 1; s <= n - no; s += off) { iv_push(&ws, (int)s); iv_push(&we, (int)(s + no)); }
        if (we.v[we.n - 1] < n) { iv_push(&ws, (int)(n - no)); iv_push(&we, n); }
    }
    int ds = 1; double rel_fs = fs;
    if (fs > 20) { ds = (int)floor(fs / 20); rel_fs = fs / ds; }
    double tol_durn = rel_fs < 10 ? 0.2 : rel_fs < 20 ? 0.1 : 0.05;
    int tol = (int)ceil(fs * tol_durn);

    ivec p = {0}, t = {0};
    double *rel = malloc(sizeof(double) * (size_t)(n / ds + 2));
    for (int w = 0; w < ws.n; w++) {
        const double *win = sig + ws.v[w] - 1;
        int W = we.v[w] - ws.v[w] + 1, Nr = 0;
        for (int i = 0; i < W; i += ds) rel[Nr++] = win[i];        /* downsample */
        msptd_scalogram(rel, Nr, rel_fs, &p, &t);
        for (int k = 0; k < p.n; k++) {
            int c = p.v[k] * ds, lo = c - tol, hi = c + tol, best;
            if (lo < 1 || hi > W) { n_bound_guards++; lo = MAXV(lo, 1); hi = MINV(hi, W); }
            best = lo; for (int i = lo + 1; i <= hi; i++) if (win[i-1] > win[best-1]) best = i;
            iv_push(peaks, best + ws.v[w] - 1);
        }
        for (int k = 0; k < t.n; k++) {
            int c = t.v[k] * ds, lo = c - tol, hi = c + tol, best;
            if (lo < 1 || hi > W) { n_bound_guards++; lo = MAXV(lo, 1); hi = MINV(hi, W); }
            best = lo; for (int i = lo + 1; i <= hi; i++) if (win[i-1] < win[best-1]) best = i;
            iv_push(onsets, best + ws.v[w] - 1);
        }
    }
    peaks->n = unique_int(peaks->v, peaks->n);
    onsets->n = unique_int(onsets->v, onsets->n);
    free(rel); free(p.v); free(t.v); free(ws.v); free(we.v);
}

/* ====================================================================== */
/* 3. detect_ppg_beats.m post-processing                                  */
/* ====================================================================== */
/* remove_repeated_peaks_and_onsets: setxor with their intersection (both sorted, unique) */
static void remove_repeated(ivec *pk, ivec *on)
{
    int i = 0, j = 0, a = 0, b = 0;
    while (i < pk->n || j < on->n) {
        if (j >= on->n || (i < pk->n && pk->v[i] < on->v[j])) pk->v[a++] = pk->v[i++];
        else if (i >= pk->n || on->v[j] < pk->v[i])          on->v[b++] = on->v[j++];
        else { i++; j++; }                                   /* in both: drop from both */
    }
    pk->n = a; on->n = b;
}

/* ensure_at_least_one_extremum_between_other_type_of_extremum
 * is_pk: other = peaks, extrema = local minima of sig; else onsets / maxima */
static void ensure_extremum_between(const double *sig, int n, ivec *other, int is_pk)
{
    int *P = calloc((size_t)n + 1, sizeof(int));          /* prefix count of strict local extrema */
    for (int i = 1; i <= n; i++) {
        int e = 0;
        if (i > 1 && i < n) e = is_pk ? (sig[i-1] < sig[i-2] && sig[i-1] < sig[i])
                                      : (sig[i-1] > sig[i-2] && sig[i-1] > sig[i]);
        P[i] = P[i-1] + e;
    }
    char *rm = malloc((size_t)other->n + 1);
    for (;;) {
        memset(rm, 0, (size_t)other->n + 1);
        int any = 0;
        for (int k = 0; k + 1 < other->n; k++) {
            int a = other->v[k], b = other->v[k+1], cnt = (b - 1 >= a + 1) ? P[b-1] - P[a] : 0;
            if (cnt == 0) {
                double va = sig[a-1], vb = sig[b-1];
                int e = is_pk ? (va <= vb ? 0 : 1) : (va >= vb ? 0 : 1);   /* min / max, first on tie */
                rm[k + e] = 1; any = 1;
            }
        }
        if (!any) break;
        int m = 0;
        for (int k = 0; k < other->n; k++) if (!rm[k]) other->v[m++] = other->v[k];
        other->n = m;
    }
    free(P); free(rm);
}

/* insert_extremum_between_other_type_of_extremum(sig, extrema, other, type)
 * is_pk: other are peaks (insert onsets at the lowest point below the line
 * joining two consecutive peaks), else other are onsets (insert peaks). */
static void insert_extremum_between(const double *sig, ivec *extrema, ivec *other, int is_pk)
{
    int tot = other->n + extrema->n;
    if (tot < 2) return;
    int *els = malloc(sizeof(int) * (size_t)tot); char *olog = malloc((size_t)tot);
    /* stable sort of [other; extrema] (values are distinct across the two sets) */
    int i = 0, j = 0, k = 0;
    while (i < other->n || j < extrema->n) {
        if (j >= extrema->n || (i < other->n && other->v[i] <= extrema->v[j])) { els[k] = other->v[i++]; olog[k++] = 1; }
        else { els[k] = extrema->v[j++]; olog[k++] = 0; }
    }
    for (int q = 0; q + 1 < tot; q++) {
        if (!(olog[q] && olog[q+1])) continue;
        int c1 = els[q], c2 = els[q+1], len = c2 - c1 + 1, temp = 1;
        double a = sig[c1-1], b = sig[c2-1], delta = (b - a) / (len - 1), best = 0;
        for (int r = 0; r < len; r++) {
            double line = (r == len - 1) ? b : a + r * delta;       /* linspace */
            double v = sig[c1 - 1 + r] - line;
            if (r == 0 || (is_pk ? v < best : v > best)) { best = v; temp = r + 1; }
        }
        if (temp == 1 || temp == len) {                    /* remove c1 from other */
            int m = 0;
            for (int z = 0; z < other->n; z++) if (other->v[z] != c1) other->v[m++] = other->v[z];
            other->n = m;
        } else {                                           /* insert into extrema, keep sorted */
            int nv = c1 - 1 + temp, pos = extrema->n;
            iv_push(extrema, nv);
            while (pos > 0 && extrema->v[pos-1] > nv) { extrema->v[pos] = extrema->v[pos-1]; pos--; }
            extrema->v[pos] = nv;
        }
    }
    free(els); free(olog);
}

static TLS int n_empty_mid = 0;             /* onset after its peak (MATLAB would error) */

static int postprocess(const double *sig, int n, ivec *pk, ivec *on, int *mid)
{
    remove_repeated(pk, on);
    ensure_extremum_between(sig, n, pk, 1);
    ensure_extremum_between(sig, n, on, 0);
    insert_extremum_between(sig, on, pk, 1);
    insert_extremum_between(sig, pk, on, 0);
    /* ensure_starts_with_onset_ends_with_peak */
    while (on->n && pk->n && on->v[0] > pk->v[0]) { memmove(pk->v, pk->v + 1, sizeof(int) * (size_t)(pk->n - 1)); pk->n--; }
    while (pk->n && on->n && pk->v[pk->n - 1] < on->v[on->n - 1]) on->n--;
    /* ensure_same_no_peaks_onsets */
    if (pk->n == 0) on->n = 0;
    if (on->n == 0) pk->n = 0;
    /* calc_mid_amp_points */
    int m = 0, nw = MINV(on->n, pk->n);
    for (int w = 0; w < nw; w++) {
        int o = on->v[w], p = pk->v[w];
        if (p < o) { n_empty_mid++; continue; }
        double des = (sig[o-1] + sig[p-1]) / 2, best = 0; int bi = o;
        for (int i = o; i <= p; i++) {
            double v = fabs(sig[i-1] - des);
            if (i == o || v < best) { best = v; bi = i; }
        }
        mid[m++] = bi;
    }
    return unique_int(mid, m);
}

/* ====================================================================== */
int msptd_nofilt(const double *ppg, int n, double fs, int *beats)
{
    ivec pk = {0}, on = {0};
    msptdfastv2(ppg, n, fs, &pk, &on);
    int m = postprocess(ppg, n, &pk, &on, beats);
    free(pk.v); free(on.v);
    return m;
}

int msptd_pb(const double *ppg, int n, double fs, int *beats)
{
    double *y = malloc(sizeof(double) * (size_t)n);
    if (!y) return -1;
    charlton_filter(ppg, n, fs, y);
    int m = msptd_nofilt(y, n, fs, beats);
    free(y);
    return m;
}

int msptd_guard_count(void) { return n_bound_guards + n_empty_mid; }
