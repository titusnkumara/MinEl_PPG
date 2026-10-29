# C implementations

**License:** `msptd_pb.c` and `msptd_pb.h` are GPL-3.0, because they translate GPL-3.0 code from the ppg-beats toolbox (`detect_ppg_beats.m`, `tidy_beats.m`; Peter H. Charlton). All other files are MIT.

Compiled ports of both detectors, for timing outside MATLAB.

| File | Port of |
| --- | --- |
| `minel_ppg.c` | `../minel_ppg.m` |
| `msptd_pb.c` | MSPTDfast v2 as run by ppg-beats: pre-filter (`preprocess_ppg_signal.m`), `msptdfastv2_beat_detector.m`, `detect_ppg_beats.m` tidying and mid-amplitude fiducial |
| `numerics.c` | MATLAB-equivalent `filtfilt`, `filter`, `median`, `prctile`, `movmedian`, `movmad` |
| `bench.c` | timing benchmark |
| `run_detector.c` | applies either detector to a binary signal file (used for verification) |
| `throughput.c` | whole-machine throughput: all recordings, all cores |
| `run_throughput.m` | packs the paper's databases and runs `throughput` from MATLAB |
| `compare_with_matlab.m` | runs the C detectors on real recordings from MATLAB and compares their beats with the MATLAB code |

## Input

`minel_ppg()` is a plain C function: `int minel_ppg(const double *ppg, int n, double fs, int *beats)`
takes the samples, their number and the sampling rate, and returns the number of beats written to
`beats` (1-based sample indices, MATLAB convention). It reads no files.

`run_detector` wraps it for use from other tools. Input file: `int32 n`, `float64 fs`, then `n`
`float64` samples; output file: `int32 m`, then `m` `int32` beat indices (little-endian).
From MATLAB:

```matlab
fid = fopen('in.bin', 'w'); fwrite(fid, numel(x), 'int32'); fwrite(fid, fs, 'double'); fwrite(fid, x, 'double'); fclose(fid);
system('./run_detector minel in.bin out.bin');
fid = fopen('out.bin'); m = fread(fid, 1, 'int32'); beats = fread(fid, m, 'int32'); fclose(fid);
```

`compare_with_matlab.m` does this for a whole database and reports any recording where the C and
MATLAB beats differ.

## Build and run

Build and run (gcc or clang, C11, POSIX threads):

```
make
./bench
```

## Verification

Both ports were checked beat for beat against the MATLAB-language reference code (run in GNU Octave 8.4 with signal and statistics packages, ppg-beats commit c397be7) on 36 synthetic recordings (64, 125 and 300 Hz; regular, AF-like and neonatal-rate rhythms; 3 min each): **36 of 36 identical for each detector.**

Details that matter for exact agreement:
- compile without fused multiply-add (`-ffp-contract=off`, set in the Makefile);
- the band-pass coefficients are computed along the same path as `butter` (prototype pole, `sftrans`, `bilinear`, `poly`), not in closed form, because a last-bit difference can change a block extremum;
- `movmad` is the *median* absolute deviation (MATLAB). Octave's `movmad` is the *mean* absolute deviation, so the Octave check used a MATLAB-equivalent `movmad`;
- the MSPTD_PB filters are designed with MATLAB's `kaiserord` + `fir1` formulas (windowed ideal sinc). Octave's `fir1` uses frequency sampling instead, so the Octave check used the C-designed taps. In MATLAB, confirm the tap counts with `length(fir1(...))`: 92 and 666 taps at 125 Hz.

**Against MATLAB R2024b on real recordings** (`compare_with_matlab.m`: MIMIC PERform Testing 125 Hz, 50 recordings; CapnoBase 300 Hz, 42; PPG-DaLiA 64 Hz, 15):
- MinEl_PPG: identical in **107 of 107** recordings;
- MSPTD_PB: identical in 104 of 107. In each of three PPG-DaLiA recordings one beat (of about 9,500-13,000) sits at a different sample, with the same beat count: a near-tie resolved differently by last-bit differences in MATLAB's `detrend`, `linspace` or `islocalmax`.

## Timing on real recordings (Apple M1 Max, one core, `compare_with_matlab.m`)

Real time factor (seconds of signal per second of computation); minimum of 5 runs per recording (3 above one hour), after one warm-up run. `compare_with_matlab.m` times the paper's workload (all recordings, each processed whole); the table below was measured on the first 50 MIMIC Testing recordings and should be replaced with a fresh run.

| Database | MinEl_PPG | MSPTD_PB (ppg-beats chain) | MSPTDfast without pre-filter |
| --- | --- | --- | --- |
| PPG-DaLiA, 64 Hz (36.0 h) | 424,000x | 29,000x (MinEl_PPG 14.5x faster) | 520,000x (MSPTDfast 1.2x faster) |
| MIMIC PERform Testing, 125 Hz (8.3 h) | 451,000x | 6,100x (MinEl_PPG 74x faster) | 544,000x (MSPTDfast 1.2x faster) |
| CapnoBase, 300 Hz (5.6 h) | 541,000x | 1,000x (MinEl_PPG 531x faster) | 376,000x (MinEl_PPG 1.4x faster) |

MinEl_PPG's cost is nearly flat across sampling rates. MSPTD_PB's grows with the sampling rate because of its long FIR pre-filter (666 taps at 125 Hz, 1,596 at 300 Hz), applied here by direct convolution as MATLAB's `filtfilt` applies it; an FFT-based filter would narrow that gap. MSPTDfast's own detector, which works on a signal decimated to about 20 Hz without filtering, costs about 20% less than MinEl_PPG at 64 and 125 Hz and more at 300 Hz.

`./bench` repeats the comparison on synthetic signals for a quick check on any machine.

## Whole-machine throughput (all cores)

`run_throughput.m` writes every recording of the selected databases to a binary pack and runs
`throughput`, which processes the recordings in parallel: each core takes the next recording
from a shared queue (longest first), and each recording is processed whole and independently,
exactly as in the paper. For each thread count (default 1, 2, 4, 8, 10) the whole set is
processed three times and the fastest wall-clock time is kept; each thread designs its filters
on a warm-up run before the clock starts. The script reports the aggregate real-time factor and
what it means in people (one person produces one second of PPG per second, so a real-time
factor of R means R people monitored in real time), the time to process one day of PPG from a
million people, and person-years processed per hour.

From a terminal: `REPS=3 ./throughput minel 1,2,4,8,10 out.bin pack1.bin [pack2.bin ...]`
(detectors: `minel`, `msptd_pb`, `msptd`). On Apple silicon the threads request the
user-initiated quality of service, which places them on performance cores first; the M1 Max
has 8 performance and 2 efficiency cores, so expect sublinear scaling from 8 to 10 threads.
The timing covers detection with signals in memory, not file reading.

### Implementation notes (all outputs bit-identical to the reference)
- the band-pass (the first-order Butterworth band-pass of the paper, `butter(1, ...)`, which as a band-pass design has two poles and three coefficients) runs in a dedicated zero-phase routine that never copies or reverses arrays; differencing and decimation share one pass, and the block extrema are found while the reverse filter pass produces its output;
- the quality gate's moving MAD uses a sliding sorted window instead of sorting each window twice;
- medians of the RR and amplitude histories are recomputed only when those histories change;
- small medians use insertion sort rather than the library qsort;
- the filter-design caches are thread-local (MSPTD_PB keeps up to eight sampling rates per thread), so both detectors can run in several threads at once with unchanged output.
