#ifndef MINEL_PPG_H
#define MINEL_PPG_H
#define MINEL_FS_MIN_DEFAULT 30.0   /* Hz: default and lowest decimated rate */
/* beats: native 1-based sample indices; buffer of at least n/2 + 1 ints.
   Returns the number of beats, or -1 on allocation failure. */
int minel_ppg(const double *ppg, int n, double fs, int *beats);
/* The same with the decimation floor fs_min (Hz, >= 30): every D-th sample
   is kept, D the largest integer with fs/D >= fs_min. Returns -2 if
   fs_min < 30. minel_ppg() is minel_ppg_fs() with fs_min = 30. */
int minel_ppg_fs(const double *ppg, int n, double fs, double fs_min, int *beats);
/* The detector's 0.5-5 Hz first-order Butterworth band-pass (butter(1, ...)):
   b and a receive its three coefficients for sampling rate fs. */
void minel_bandpass(double fs, double b[3], double a[3]);
#endif
