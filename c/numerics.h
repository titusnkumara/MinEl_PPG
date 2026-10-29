#ifndef NUMERICS_H
#define NUMERICS_H
void   sort_dbl(double *x, int n);
void   sort_int(int *x, int n);
int    unique_int(int *x, int n);
double median_dbl(const double *x, int n);
double prctile_dbl(const double *x, int n, double p);
void   movmedian_odd(const double *x, int n, int w, double *out);
void   movmad_even(const double *x, int n, int w, double *out);
int    filtfilt(const double *b, const double *a, int nc, const double *x, int n, double *out);
int filtfilt3_extrema(const double b[3], const double a[3], const double *x, int n, double *out, const long *blk_e, int nb, int *imin, int *imax);
#endif
