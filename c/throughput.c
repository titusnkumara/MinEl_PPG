/* throughput: whole-machine throughput of a beat detector, using all cores.
 *
 *   throughput [-f FS_MIN] DET THREADS out.bin pack1.bin [pack2.bin ...]
 *
 *   -f FS_MIN MinEl_PPG's decimation floor in Hz (default 30, minimum 30)
 *   DET      minel | msptd_pb | msptd | msptd_bp
 *            msptd    = MSPTDfast without pre-filter
 *            msptd_bp = MSPTDfast after the same band-pass MinEl_PPG uses
 *                       (0.5-5 Hz first-order Butterworth, zero-phase, at the native rate)
 *   THREADS  comma-separated thread counts, e.g. 1,2,4,8,10 ("all" = every core)
 *   out.bin  results for MATLAB: int32 k, float64 signal_s, then k x (float64 threads, float64 wall_s)
 *   pack     int32 nrec, then per recording: float64 fs, int32 n, float64 x[n]
 *            (write it with run_throughput.m; recordings may have different fs)
 *
 * Every recording is processed whole and independently, exactly as in the
 * paper; recordings are handed to the threads from a shared queue, longest
 * first, so no core idles while long recordings remain. For each thread
 * count the whole set is processed REPS times (default 3, environment
 * variable REPS) and the fastest wall-clock time is kept. Each worker fills
 * its own filter-design caches on a short warm-up before the clock starts,
 * so design is excluded for every arm, as in the paper. Signals are held in
 * memory: the timing covers detection only, not file reading.
 *
 * Throughput is reported as an aggregate real-time factor: seconds of PPG
 * processed per second of wall-clock time. Because one person produces one
 * second of PPG per second, this is also the number of people whose PPG
 * the machine could keep up with in real time (offline, whole-recording
 * processing; no I/O or streaming latency included).                      */
#define _POSIX_C_SOURCE 200809L
#define _DARWIN_C_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <pthread.h>
#include <stdatomic.h>
#include <unistd.h>
#ifdef __APPLE__
#include <pthread/qos.h>
#endif
#include "minel_ppg.h"
#include "msptd_pb.h"
#include "numerics.h"

typedef int (*det_fn)(const double *, int, double, int *);
typedef struct { double fs; int n; double *x; } rec_t;
static double g_fs_min = MINEL_FS_MIN_DEFAULT;
static int minel_det(const double *x, int n, double fs, int *b) { return minel_ppg_fs(x, n, fs, g_fs_min, b); }

/* MSPTDfast given the simplest band-pass: MinEl_PPG's own filter, applied
   to the raw PPG at the native rate, then MSPTDfast without its FIR pre-filter */
static int msptd_bp(const double *x, int n, double fs, int *beats)
{
    double b[3], a[3];
    minel_bandpass(fs, b, a);
    double *y = malloc(sizeof(double) * (size_t)n);
    if (!y) return -1;
    if (filtfilt3_extrema(b, a, x, n, y, NULL, 0, NULL, NULL) != 0) memcpy(y, x, sizeof(double) * (size_t)n);
    int m = msptd_nofilt(y, n, fs, beats);
    free(y);
    return m;
}

static double now_s(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return t.tv_sec + 1e-9 * t.tv_nsec; }

/* ---------------------------------------------------------------- data */
static rec_t *R = NULL; static int NR = 0, NMAX = 0;
static int load_pack(const char *path)
{
    FILE *f = fopen(path, "rb"); if (!f) { fprintf(stderr, "cannot open %s\n", path); return -1; }
    int m; if (fread(&m, 4, 1, f) != 1) { fclose(f); return -1; }
    R = realloc(R, sizeof(rec_t) * (size_t)(NR + m));
    for (int i = 0; i < m; i++) {
        rec_t *r = &R[NR];
        if (fread(&r->fs, 8, 1, f) != 1 || fread(&r->n, 4, 1, f) != 1) { fclose(f); return -1; }
        r->x = malloc(sizeof(double) * (size_t)r->n);
        if (!r->x || fread(r->x, 8, (size_t)r->n, f) != (size_t)r->n) { fclose(f); return -1; }
        if (r->n > NMAX) NMAX = r->n;
        NR++;
    }
    fclose(f);
    return m;
}
static int cmp_len_desc(const void *a, const void *b)
{
    const rec_t *x = a, *y = b;
    double cx = (double)x->n, cy = (double)y->n;          /* longest first */
    return (cx < cy) - (cx > cy);
}

/* ------------------------------------------------------------- workers */
typedef struct {
    det_fn f;
    atomic_int next;
    pthread_mutex_t mu; pthread_cond_t cv;
    int ready, go;
    double t_end_max;
} run_t;

static void *worker(void *arg)
{
    run_t *J = arg;
#ifdef __APPLE__
    pthread_set_qos_class_self_np(QOS_CLASS_USER_INITIATED, 0);   /* prefer performance cores */
#endif
    int *beats = malloc(sizeof(int) * (size_t)(NMAX + 1));
    /* warm-up: one short run per distinct sampling rate fills this thread's caches */
    for (int i = 0; i < NR; i++) {
        int seen = 0;
        for (int j = 0; j < i; j++) if (R[j].fs == R[i].fs) { seen = 1; break; }
        if (seen) continue;
        int nw = R[i].n < (int)(120 * R[i].fs) ? R[i].n : (int)(120 * R[i].fs);
        J->f(R[i].x, nw, R[i].fs, beats);
    }
    pthread_mutex_lock(&J->mu);
    J->ready++; pthread_cond_broadcast(&J->cv);
    while (!J->go) pthread_cond_wait(&J->cv, &J->mu);
    pthread_mutex_unlock(&J->mu);

    for (;;) {
        int i = atomic_fetch_add(&J->next, 1);
        if (i >= NR) break;
        J->f(R[i].x, R[i].n, R[i].fs, beats);
    }
    double t_end = now_s();
    pthread_mutex_lock(&J->mu);
    if (t_end > J->t_end_max) J->t_end_max = t_end;
    pthread_mutex_unlock(&J->mu);
    free(beats);
    msptd_release_cache();
    return NULL;
}

/* wall-clock time to process every recording once with T threads */
static double run_once(det_fn f, int T)
{
    run_t J; memset(&J, 0, sizeof J);
    J.f = f; atomic_init(&J.next, 0);
    pthread_mutex_init(&J.mu, NULL); pthread_cond_init(&J.cv, NULL);
    pthread_t *th = malloc(sizeof(pthread_t) * (size_t)T);
    for (int t = 0; t < T; t++) pthread_create(&th[t], NULL, worker, &J);
    pthread_mutex_lock(&J.mu);
    while (J.ready < T) pthread_cond_wait(&J.cv, &J.mu);
    double t0 = now_s();
    J.go = 1; pthread_cond_broadcast(&J.cv);
    pthread_mutex_unlock(&J.mu);
    for (int t = 0; t < T; t++) pthread_join(th[t], NULL);
    free(th);
    pthread_mutex_destroy(&J.mu); pthread_cond_destroy(&J.cv);
    return J.t_end_max - t0;
}

/* ----------------------------------------------------------------- main */
int main(int argc, char **argv)
{
    if (argc >= 3 && !strcmp(argv[1], "-f")) {              /* optional decimation floor */
        g_fs_min = atof(argv[2]);
        if (!(g_fs_min >= MINEL_FS_MIN_DEFAULT)) { fprintf(stderr, "-f: FS_MIN must be >= %g Hz\n", MINEL_FS_MIN_DEFAULT); return 1; }
        argv += 2; argc -= 2;
    }
    if (argc < 5) {
        fprintf(stderr, "usage: throughput [-f FS_MIN] minel|msptd_pb|msptd|msptd_bp THREADS out.bin pack1.bin [pack2.bin ...]\n");
        return 1;
    }
    det_fn f = !strcmp(argv[1], "minel") ? minel_det : !strcmp(argv[1], "msptd_pb") ? msptd_pb :
               !strcmp(argv[1], "msptd") ? msptd_nofilt : !strcmp(argv[1], "msptd_bp") ? msptd_bp : NULL;
    if (!f) { fprintf(stderr, "unknown detector %s\n", argv[1]); return 1; }
    for (int a = 4; a < argc; a++) if (load_pack(argv[a]) < 0) { fprintf(stderr, "bad pack %s\n", argv[a]); return 2; }
    if (NR == 0) return 2;
    qsort(R, (size_t)NR, sizeof(rec_t), cmp_len_desc);
    double S = 0; for (int i = 0; i < NR; i++) S += R[i].n / R[i].fs;

    int ncpu = (int)sysconf(_SC_NPROCESSORS_ONLN), T[64], nt = 0;
    if (!strcmp(argv[2], "all")) T[nt++] = ncpu;
    else for (char *p = strtok(argv[2], ","); p && nt < 64; p = strtok(NULL, ",")) T[nt++] = atoi(p);
    int reps = getenv("REPS") ? atoi(getenv("REPS")) : 3; if (reps < 1) reps = 1;

    printf("%s: %d recordings, %.1f h of PPG, %d logical cores, best of %d runs",
           argv[1], NR, S / 3600, ncpu, reps);
    if (f == minel_det) printf(", decimation floor %g Hz", g_fs_min);
    printf("\n");
    printf("%8s %12s %16s %9s %11s\n", "threads", "wall (s)", "real-time factor", "speed-up", "efficiency");
    double *W = malloc(sizeof(double) * (size_t)nt), w1 = 0;
    for (int k = 0; k < nt; k++) {
        double best = 1e30;
        for (int r = 0; r < reps; r++) { double w = run_once(f, T[k]); if (w < best) best = w; }
        W[k] = best;
        if (T[k] == 1) w1 = best;
        if (w1 > 0) printf("%8d %12.4f %16.0f %8.2fx %10.0f%%\n", T[k], best, S / best, w1 / best, 100 * w1 / best / T[k]);
        else        printf("%8d %12.4f %16.0f\n", T[k], best, S / best);
    }
    FILE *o = fopen(argv[3], "wb");
    if (o) {
        fwrite(&nt, 4, 1, o); fwrite(&S, 8, 1, o);
        for (int k = 0; k < nt; k++) { double tk = T[k]; fwrite(&tk, 8, 1, o); fwrite(&W[k], 8, 1, o); }
        fclose(o);
    }
    return 0;
}
