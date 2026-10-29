/* bench: time MinEl_PPG against MSPTDfast (ppg-beats chain) in C.
 * Synthetic PPG (pulse + dicrotic wave + noise + baseline wander), 72 bpm
 * with beat-to-beat variability. Filter design is excluded by a warm-up run
 * (cached, as in the paper). Reported: minimum of REPS runs, single thread.
 *   bench [-f FS_MIN]   MinEl_PPG's decimation floor in Hz (default 30, minimum 30) */
#define _POSIX_C_SOURCE 199309L
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846   /* not defined in strict C/POSIX modes (e.g. macOS) */
#endif
#include <time.h>
#include <string.h>
#include "minel_ppg.h"
#include "msptd_pb.h"

static double now(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return t.tv_sec + 1e-9 * t.tv_nsec; }
static double randn(void)
{
    double u1 = (rand() + 1.0) / (RAND_MAX + 2.0), u2 = (rand() + 1.0) / (RAND_MAX + 2.0);
    return sqrt(-2 * log(u1)) * cos(2 * M_PI * u2);
}
static double *synth(double fs, double dur, int *n, int *nbeats)
{
    *n = (int)(fs * dur) + 1;
    double *x = calloc((size_t)*n, sizeof(double)), t = 0.5;
    int nb = 0;
    while (t < dur - 0.5) {
        double rr = 0.83 + 0.04 * randn(), a = 1 + 0.2 * randn();
        int c = (int)(t * fs), w = (int)(0.6 * fs);
        for (int i = c - w; i <= c + w; i++) if (i >= 0 && i < *n) {
            double s = i / fs - t;
            x[i] += a * exp(-pow(s / 0.09, 2)) + 0.3 * a * exp(-pow((s - 0.3) / 0.07, 2));
        }
        t += rr; nb++;
    }
    for (int i = 0; i < *n; i++) x[i] += 0.03 * randn() + 0.2 * sin(2 * M_PI * 0.2 * i / fs);
    *nbeats = nb;
    return x;
}
typedef int (*det_fn)(const double *, int, double, int *);
static double g_fs_min = MINEL_FS_MIN_DEFAULT;
static int minel_det(const double *x, int n, double fs, int *b) { return minel_ppg_fs(x, n, fs, g_fs_min, b); }
static double time_it(det_fn f, const double *x, int n, double fs, int reps, int *nout)
{
    int *b = malloc(sizeof(int) * (size_t)(n + 1));
    *nout = f(x, n, fs, b);                       /* warm-up: caches filter design */
    double best = 1e30;
    for (int r = 0; r < reps; r++) { double t0 = now(); f(x, n, fs, b); double dt = now() - t0; if (dt < best) best = dt; }
    free(b);
    return best;
}
int main(int argc, char **argv)
{
    if (argc >= 3 && !strcmp(argv[1], "-f")) {
        g_fs_min = atof(argv[2]);
        if (!(g_fs_min >= MINEL_FS_MIN_DEFAULT)) { fprintf(stderr, "-f: FS_MIN must be >= %g Hz\n", MINEL_FS_MIN_DEFAULT); return 1; }
    }
    printf("MinEl_PPG decimation floor: %g Hz\n", g_fs_min);
    srand(2024);
    double fss[] = {64, 125, 300}, durs[] = {600, 3600};
    printf("%-6s %-7s | %11s %11s %11s | %9s %9s | %s\n", "fs", "length", "MinEl_PPG", "MSPTD_PB", "MSPTD", "vs PB", "vs MSPTD", "beats (true) U / PB / MSPTD");
    for (int d = 0; d < 2; d++) for (int k = 0; k < 3; k++) {
        int n, nb, nu, np, nm;
        double fs = fss[k], *x = synth(fs, durs[d], &n, &nb);
        int reps = durs[d] > 1000 ? 3 : 7;
        double tu = time_it(minel_det, x, n, fs, reps, &nu);
        double tp = time_it(msptd_pb, x, n, fs, reps, &np);
        double tm = time_it(msptd_nofilt, x, n, fs, reps, &nm);
        printf("%-3.0f Hz %5.0f s | %9.2f ms %8.1f ms %8.1f ms | %8.1fx %8.1fx | %d (%d) / %d / %d\n",
               fs, durs[d], 1e3 * tu, 1e3 * tp, 1e3 * tm, tp / tu, tm / tu, nu, nb, np, nm);
        printf("%-15s | RTF %8.0fk %8.0fk %8.0fk\n", "", durs[d] / tu / 1e3, durs[d] / tp / 1e3, durs[d] / tm / 1e3);
        free(x);
    }
    return 0;
}
