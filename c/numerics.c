/* numerics.c - MATLAB-equivalent numerics shared by both detectors.
 *
 * Operation order follows MATLAB/Octave so that results match to the last
 * bit where possible: filter() is direct-form II transposed, filtfilt() uses
 * 3*(order) samples of odd reflection at each end and steady-state initial
 * conditions, median/prctile/movmedian/movmad follow MATLAB's definitions.
 * Compile with -ffp-contract=off (no fused multiply-add).
 */
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "numerics.h"

static int cmp_dbl(const void *a, const void *b)
{
    double x = *(const double *)a, y = *(const double *)b;
    return (x > y) - (x < y);
}
static int cmp_int(const void *a, const void *b)
{
    int x = *(const int *)a, y = *(const int *)b;
    return (x > y) - (x < y);
}

void sort_dbl(double *x, int n) { qsort(x, (size_t)n, sizeof(double), cmp_dbl); }
void sort_int(int *x, int n)    { qsort(x, (size_t)n, sizeof(int), cmp_int); }

int unique_int(int *x, int n)          /* sort + unique in place, returns new length */
{
    if (n <= 0) return 0;
    sort_int(x, n);
    int m = 1;
    for (int i = 1; i < n; i++) if (x[i] != x[m-1]) x[m++] = x[i];
    return m;
}

double median_dbl(const double *x, int n)   /* MATLAB median: mean of middle two if even */
{
    if (n <= 0) return NAN;
    double buf[64], *s = n <= 64 ? buf : malloc((size_t)n * sizeof(double));
    memcpy(s, x, (size_t)n * sizeof(double));
    if (n <= 64) {                                    /* insertion sort: small windows */
        for (int i = 1; i < n; i++) {
            double v = s[i]; int j = i - 1;
            while (j >= 0 && s[j] > v) { s[j+1] = s[j]; j--; }
            s[j+1] = v;
        }
    } else sort_dbl(s, n);
    double m = (n % 2) ? s[n/2] : (s[n/2 - 1] + s[n/2]) / 2.0;
    if (s != buf) free(s);
    return m;
}

/* MATLAB prctile (default method): sorted values sit at percentiles
 * 100*(i-0.5)/n; linear interpolation between them, clamped at the ends. */
double prctile_dbl(const double *x, int n, double p)
{
    if (n <= 0) return NAN;
    double *s = malloc((size_t)n * sizeof(double));
    memcpy(s, x, (size_t)n * sizeof(double));
    sort_dbl(s, n);
    double pos = n * p / 100.0 + 0.5, r;
    if (pos <= 1.0) r = s[0];
    else if (pos >= n) r = s[n-1];
    else {
        int i = (int)floor(pos);
        double f = pos - i;
        r = s[i-1] + f * (s[i] - s[i-1]);
    }
    free(s);
    return r;
}

/* movmedian with odd window w, endpoints shrink */
void movmedian_odd(const double *x, int n, int w, double *out)
{
    int h = w / 2;
    for (int i = 0; i < n; i++) {
        int a = i - h < 0 ? 0 : i - h, b = i + h > n - 1 ? n - 1 : i + h;
        out[i] = median_dbl(x + a, b - a + 1);
    }
}

/* movmad(x, w) with even w: window [i - w/2, i + w/2 - 1], endpoints shrink;
 * mad = median(|x - median(x)|) (not scaled) */
void movmad_even(const double *x, int n, int w, double *out)
{
    int kb = w / 2, kf = w / 2 - 1;
    double tmp[64];
    for (int i = 0; i < n; i++) {
        int a = i - kb < 0 ? 0 : i - kb, b = i + kf > n - 1 ? n - 1 : i + kf;
        int m = b - a + 1;
        double med = median_dbl(x + a, m);
        for (int j = 0; j < m; j++) tmp[j] = fabs(x[a+j] - med);
        out[i] = median_dbl(tmp, m);
    }
}

/* filter(b, a, x, zi) - direct form II transposed, a[0] == 1, exactly as
 * liboctave: y = z0 + b0*x;  z_j = (z_{j+1} - a_{j+1}*y) + b_{j+1}*x;
 * z_last = b_n*x - a_n*y. zi (length nc-1) is overwritten. */
static void filter_df2t(const double *b, const double *a, int nc,
                        const double *x, int n, double *y, double *z)
{
    int ns = nc - 1;
    for (int i = 0; i < n; i++) {
        double xi = x[i];
        double yi = z[0] + b[0] * xi;
        for (int j = 0; j < ns - 1; j++)
            z[j] = z[j+1] - a[j+1] * yi + b[j+1] * xi;
        z[ns-1] = b[ns] * xi - a[ns] * yi;
        y[i] = yi;
    }
}

/* FIR-only filter (a = [1 0 ... 0]) with the same arithmetic order as the
 * general DF2T path but without the O(order) state update per sample: the
 * output is built innermost-first, y[t] = (((z_{t+1} + b_t x0) + ...) +
 * b_1 x_{t-1}) + b_0 x_t, which is the DF2T recursion unrolled. */
static void filter_fir(const double *b, int nc, const double *x, int n,
                       double *y, const double *zi)
{
    int ns = nc - 1;
    for (int t = 0; t < n; t++) {
        int jmax = t < ns ? t : ns;
        double acc;
        if (t < ns) acc = zi[t];                 /* z_{t+1}(initial) */
        else        acc = b[ns] * x[t - ns] - 0.0;
        int jstart = (t < ns) ? jmax : ns - 1;
        for (int j = jstart; j >= 1; j--)
            acc = acc - 0.0 + b[j] * x[t - j];
        y[t] = acc + b[0] * x[t];
    }
}

/* filtfilt(b, a, x) as Octave/MATLAB: odd reflection of 3*(nc-1) samples at
 * each end, steady-state initial conditions zi_m = sum_{j>=m} (b_j - kdc a_j),
 * forward pass with zi*v(1), reverse pass with zi*v(end). a == NULL means FIR. */
int filtfilt(const double *b, const double *a_in, int nc, const double *x, int n, double *out)
{
    int lr = 3 * (nc - 1);
    if (n <= lr) return -1;
    double *a = calloc((size_t)nc, sizeof(double));
    if (a_in) memcpy(a, a_in, (size_t)nc * sizeof(double)); else a[0] = 1.0;

    double sb = 0, sa = 0;
    for (int j = 0; j < nc; j++) { sb += b[j]; sa += a[j]; }
    double kdc = sb / sa;
    double *si = malloc((size_t)nc * sizeof(double));
    double c = 0;                                /* cumsum from the end */
    for (int j = nc - 1; j >= 0; j--) { c = c + (b[j] - kdc * a[j]); si[j] = c; }
    /* si(1) removed: zi[m-1] = si[m], m = 1..nc-1 */

    int nv = n + 2 * lr;
    double *v = malloc((size_t)nv * sizeof(double)), *w = malloc((size_t)nv * sizeof(double));
    double *z = malloc((size_t)nc * sizeof(double));
    for (int i = 0; i < lr; i++) v[i] = 2 * x[0] - x[lr - i];
    memcpy(v + lr, x, (size_t)n * sizeof(double));
    for (int i = 0; i < lr; i++) v[lr + n + i] = 2 * x[n-1] - x[n - 2 - i];

    for (int m = 0; m < nc - 1; m++) z[m] = si[m+1] * v[0];
    if (a_in) filter_df2t(b, a, nc, v, nv, w, z); else filter_fir(b, nc, v, nv, w, z);
    for (int i = 0; i < nv / 2; i++) { double t = w[i]; w[i] = w[nv-1-i]; w[nv-1-i] = t; }
    for (int m = 0; m < nc - 1; m++) z[m] = si[m+1] * w[0];
    if (a_in) filter_df2t(b, a, nc, w, nv, v, z); else filter_fir(b, nc, w, nv, v, z);
    for (int i = 0; i < n; i++) out[i] = v[nv - 1 - lr - i];

    free(a); free(si); free(v); free(w); free(z);
    return 0;
}

/* filtfilt for a 3-coefficient IIR (nc = 3) - e.g. the first-order Butterworth
 * band-pass butter(1, [f1 f2], 'bandpass'), which has two poles - with the arithmetic of filtfilt()
 * above, but without copying or reversing arrays: the forward pass writes w,
 * the reverse pass walks w backwards and writes out[] directly. If blk_e is
 * given, the minimum and maximum of each block (samples blk_e[q]+1 ..
 * blk_e[q+1], 1-based) are found while the reverse pass produces them;
 * scanning backwards with <= / >= keeps the earliest sample on ties, as
 * MATLAB's min/max do. */
int filtfilt3_extrema(const double b[3], const double a[3], const double *x, int n, double *out,
                      const long *blk_e, int nb, int *imin, int *imax)
{
    const int lr = 6;                               /* 3 * (nc - 1) */
    if (n <= lr) return -1;
    const double kdc = (b[0] + b[1] + b[2]) / (a[0] + a[1] + a[2]);
    const double c2 = 0.0 + (b[2] - kdc * a[2]);     /* cumsum from the end */
    const double c1 = c2 + (b[1] - kdc * a[1]);      /* zi = [c1 c2] */
    const int nv = n + 2 * lr;
    double *w = malloc(sizeof(double) * (size_t)nv);
    if (!w) return -2;
    const double b0 = b[0], b1 = b[1], b2 = b[2], a1 = a[1], a2 = a[2];

    /* forward pass over v = [reflection; x; reflection], generated on the fly */
    double v0 = 2 * x[0] - x[lr], z0 = c1 * v0, z1 = c2 * v0, yi, xi;
    for (int i = 0; i < lr; i++) {
        xi = 2 * x[0] - x[lr - i];
        yi = z0 + b0 * xi; z0 = z1 - a1 * yi + b1 * xi; z1 = b2 * xi - a2 * yi; w[i] = yi;
    }
    for (int i = 0; i < n; i++) {
        xi = x[i];
        yi = z0 + b0 * xi; z0 = z1 - a1 * yi + b1 * xi; z1 = b2 * xi - a2 * yi; w[lr + i] = yi;
    }
    for (int i = 0; i < lr; i++) {
        xi = 2 * x[n-1] - x[n - 2 - i];
        yi = z0 + b0 * xi; z0 = z1 - a1 * yi + b1 * xi; z1 = b2 * xi - a2 * yi; w[lr + n + i] = yi;
    }

    /* reverse pass: trailing pad, signal (with block extrema), leading pad */
    double wl = w[nv-1];
    z0 = c1 * wl; z1 = c2 * wl;
    for (int i = nv - 1; i >= lr + n; i--) {
        xi = w[i]; yi = z0 + b0 * xi; z0 = z1 - a1 * yi + b1 * xi; z1 = b2 * xi - a2 * yi;
    }
    long last = blk_e ? blk_e[nb] : 0;               /* samples beyond the last block are not coded */
    int q = nb - 1; double vmin = 0, vmax = 0;
    for (int j = n - 1; j >= 0; j--) {
        xi = w[lr + j]; yi = z0 + b0 * xi; z0 = z1 - a1 * yi + b1 * xi; z1 = b2 * xi - a2 * yi;
        out[j] = yi;
        long s1 = j + 1;
        if (blk_e && s1 <= last) {
            while (s1 <= blk_e[q]) q--;
            if (s1 == blk_e[q+1]) { imin[q] = imax[q] = (int)s1; vmin = vmax = yi; }
            else {
                if (yi <= vmin) { vmin = yi; imin[q] = (int)s1; }
                if (yi >= vmax) { vmax = yi; imax[q] = (int)s1; }
            }
        }
    }
    free(w);
    return 0;
}
