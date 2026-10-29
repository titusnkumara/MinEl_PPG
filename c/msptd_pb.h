/* GPL-3.0: see msptd_pb.c for license and attribution. */
#ifndef MSPTD_PB_H
#define MSPTD_PB_H
/* MSPTDfast v2 as run by ppg-beats: pre-filter, detector, onset-peak tidying,
   mid-amplitude fiducial. beats: 1-based indices, buffer >= n ints. */
int msptd_pb(const double *ppg, int n, double fs, int *beats);
/* the same without the pre-filter */
int msptd_nofilt(const double *ppg, int n, double fs, int *beats);
/* filter taps used for fs (design cached); returns tap counts */
void msptd_filter_taps(double fs, const double **b_high, int *n_high, const double **b_low, int *n_low);

/* boundary guards hit by the calling thread */
int msptd_guard_count(void);
/* free the calling thread's filter-design cache (call before a worker thread exits) */
void msptd_release_cache(void);
#endif
