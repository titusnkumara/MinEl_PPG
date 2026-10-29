/* run_detector: apply a detector to a binary signal file.
 *   run_detector [-f FS_MIN] minel|msptd_pb  in.bin  out.bin   (beats)
 *   run_detector [-f FS_MIN] time in.bin [out.bin]              (timing: s per detector; 3 x float64 to out.bin)
 * -f FS_MIN : MinEl_PPG's decimation floor in Hz (default 30, minimum 30)
 * in.bin : int32 n, float64 fs, float64 x[n]
 * out.bin: int32 m, int32 beats[m] (1-based sample indices)            */
#define _POSIX_C_SOURCE 199309L   /* clock_gettime on macOS and Linux */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "minel_ppg.h"
#include "msptd_pb.h"
#include <time.h>
static double now_s(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return t.tv_sec + 1e-9 * t.tv_nsec; }
typedef int (*det_fn)(const double *, int, double, int *);
static double g_fs_min = MINEL_FS_MIN_DEFAULT;
static int minel_det(const double *x, int n, double fs, int *b) { return minel_ppg_fs(x, n, fs, g_fs_min, b); }
/* minimum of reps runs after one warm-up run (fills the filter-design caches) */
static double time_min(det_fn f, const double *x, int n, double fs, int *buf, int reps)
{
    f(x, n, fs, buf);
    double best = 1e30;
    for (int r = 0; r < reps; r++) { double t0 = now_s(); f(x, n, fs, buf); double dt = now_s() - t0; if (dt < best) best = dt; }
    return best;
}
int main(int argc, char **argv)
{
    if (argc >= 3 && !strcmp(argv[1], "-f")) {              /* optional decimation floor */
        g_fs_min = atof(argv[2]);
        if (!(g_fs_min >= MINEL_FS_MIN_DEFAULT)) { fprintf(stderr, "-f: FS_MIN must be >= %g Hz\n", MINEL_FS_MIN_DEFAULT); return 1; }
        argv += 2; argc -= 2;
    }
    if (argc == 4 && !strcmp(argv[1], "taps")) {        /* run_detector taps fs out.bin */
        const double *bh, *bl; int nh, nl; double fs = atof(argv[2]);
        msptd_filter_taps(fs, &bh, &nh, &bl, &nl);
        FILE *g = fopen(argv[3], "wb");
        fwrite(&nh, 4, 1, g); fwrite(bh, 8, (size_t)nh, g); fwrite(&nl, 4, 1, g); fwrite(bl, 8, (size_t)nl, g);
        fclose(g); return 0;
    }
    if ((argc == 3 || argc == 4) && !strcmp(argv[1], "time")) {   /* run_detector time in.bin [out.bin] */
        FILE *g = fopen(argv[2], "rb"); if (!g) return 2;
        int n; double fs;
        if (fread(&n, 4, 1, g) != 1 || fread(&fs, 8, 1, g) != 1) return 3;
        double *x = malloc(sizeof(double) * (size_t)n);
        if (fread(x, 8, (size_t)n, g) != (size_t)n) return 3;
        fclose(g);
        int *b = malloc(sizeof(int) * (size_t)(n + 1));
        int reps = n / fs > 3600 ? 3 : 5;
        /* seconds: MinEl_PPG, MSPTD_PB, MSPTDfast without pre-filter */
        double t[3] = { time_min(minel_det, x, n, fs, b, reps),
                        time_min(msptd_pb, x, n, fs, b, reps),
                        time_min(msptd_nofilt, x, n, fs, b, reps) };
        if (argc == 4) {                                   /* 3 x float64 to file */
            FILE *o = fopen(argv[3], "wb"); if (!o) return 4;
            fwrite(t, 8, 3, o); fclose(o);
        } else {
            printf("%.9f %.9f %.9f\n", t[0], t[1], t[2]);
        }
        free(x); free(b); return 0;
    }
    if (argc != 4) { fprintf(stderr, "usage: run_detector [-f FS_MIN] minel|msptd_pb in.bin out.bin\n"); return 1; }
    FILE *f = fopen(argv[2], "rb"); if (!f) return 2;
    int n; double fs;
    if (fread(&n, 4, 1, f) != 1 || fread(&fs, 8, 1, f) != 1) return 3;
    double *x = malloc(sizeof(double) * (size_t)n);
    if (fread(x, 8, (size_t)n, f) != (size_t)n) return 3;
    fclose(f);
    int *b = malloc(sizeof(int) * (size_t)(n + 1)), m;
    if (!strcmp(argv[1], "minel")) m = minel_det(x, n, fs, b);
    else m = msptd_pb(x, n, fs, b);
    if (msptd_guard_count()) fprintf(stderr, "note: %d boundary guards\n", msptd_guard_count());
    f = fopen(argv[3], "wb"); fwrite(&m, 4, 1, f); fwrite(b, 4, (size_t)(m > 0 ? m : 0), f); fclose(f);
    return 0;
}
