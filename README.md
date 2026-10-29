# MinEl_PPG

Fast, accurate photoplethysmogram (PPG) beat detection from **sparse extremum coding**. The input is reduced to about 30 Hz, the band-passed PPG derivative is split into 0.2 s blocks, and only each block's minimum and maximum are kept, with their sample indices. A rhythm-aware detector works on this stream of ten values per second, at any sampling rate, and each beat's timing is refined back to native resolution.

This repository contains the detector and the full code to reproduce the paper:

> **MinEl_PPG: Sparse Min–Max Coding for Fast, Accurate PPG Beat Detection.** [authors, journal, year, DOI ____]

## Use the detector

`minel_ppg.m` is self-contained (one file). Copy it into your project:

```matlab
beats = minel_ppg(ppg, fs);   % raw PPG and sampling rate (Hz) in, beat sample indices out
```

Each beat is reported at the point of steepest upslope, refined by parabolic interpolation and mapped to the nearest native sample. Requires MATLAB with the Signal Processing Toolbox (`butter`, `filtfilt`) and the Statistics and Machine Learning Toolbox (`prctile`). See `example.m`.

The input is reduced by keeping every D-th sample, where D is the largest integer that leaves at least 30 Hz (64 → 32 Hz, 125 → 31.25 Hz, 300 → 30 Hz). To keep more samples, raise the floor:

```matlab
beats = minel_ppg(ppg, fs, struct('fs_min', 100));   % floor in Hz, at least 30
```

A C version with the same output is in `c/`: `minel_ppg(ppg, n, fs, beats)` uses the 30 Hz default, and `minel_ppg_fs(ppg, n, fs, fs_min, beats)` takes the floor as an argument.

## Reproduce the paper

```
MinEl_PPG/
├── minel_ppg.m           the detector
├── example.m             minimal usage example
├── download_data.m       fetches the seven databases from Zenodo
├── run_paper.m           reproduces every table and figure number
├── benchmark/
│   ├── load_databases.m  loaders for the seven databases
│   ├── patient_overlap.m held-out mask: removes patients shared with other MIMIC PERform datasets
│   ├── ecg_references.m  ECG reference beats (Pan–Tompkins, MinEl, TiTen)
│   ├── pan_tompkin.m     Pan–Tompkins implementation (third party, see header)
│   ├── comparators.m     MSPTDfast with the ppg-beats processing chain (pre-filter, tidying, mid-amplitude fiducial)
│   ├── metrics.m         matching, alignment, F1, heart rate and HRV
│   └── run_statistics.m  every statistical analysis in the paper
└── c/                    C ports of both detectors, verified against the MATLAB code, and the
                          whole-machine throughput benchmark (see c/README.md)
```

1. **Requirements.** MATLAB R2024b, Signal Processing Toolbox, Statistics and Machine Learning Toolbox, and the [ppg-beats toolbox](https://github.com/peterhcharlton/ppg-beats) version ____ on the MATLAB path (it provides MSPTDfast). For the C benchmarks: gcc or clang (C11) and POSIX threads.
2. **Data.** Run `download_data`. It fetches the seven databases from Zenodo ([10.5281/zenodo.23721153](https://doi.org/10.5281/zenodo.23721153), 437 MB), checks each file's MD5, resumes interrupted downloads and unzips PPG-DaLiA. By default it writes to the parent of the repository folder, which is where `run_paper.m` looks (`DATA = '../'`):
   ```
   ../mimic-large-test/   ../mimic_test/   ../mimic_train/   ../mimic_perform_af/
   ../bidmc/              ../capnobase/    ../dalia/dalia_mat/
   ```
   Use `download_data(folder)` and set `DATA` in `run_paper.m` to put them elsewhere. The databases are redistributed from their original sources, whose licences and citations apply: MIMIC PERform Large Testing, Testing, Training and AF (Charlton et al 2022, via ppg-beats), BIDMC (Pimentel et al 2017), CapnoBase (Karlen et al 2013) and PPG-DaLiA (Reiss et al 2019).
3. **Held-out set.** Nothing to do: whenever Large Testing is loaded, records from patients who appear in any other MIMIC PERform dataset are removed automatically (matched on MIMIC-III subject ID), and the overlap counts are printed (796 recordings from 111 patients remain).
4. **Run** `run_paper`. Every detector processes each recording whole, identically on all seven databases. It prints the per-database tables and every statistical analysis, and saves all results to `ppg_benchmark_ACC.mat`. Random seeds are fixed.
5. **Throughput (optional).** In `c/`, run `make`, then `run_throughput` in MATLAB. It processes all recordings on every core and reports the aggregate real-time factor and the number of machines needed to analyse PPG from a given population in real time. Its `FS_MIN` setting must match the detector's floor (30 Hz, as in the paper).

Timings depend on the machine; the paper's were measured in MATLAB and, for whole-machine throughput, in C, on an Apple M1 Max (10 cores: 8 performance, 2 efficiency; 32 GB RAM).

## Citation

If you use this code, please cite the paper above. If you use the data, also cite the Zenodo record and the original databases.

## License

MIT (see `LICENSE`), with two exceptions: `benchmark/pan_tompkin.m` keeps its original license (see its header), and `c/msptd_pb.c` / `c/msptd_pb.h` are GPL-3.0 because they translate GPL-3.0 code from the ppg-beats toolbox; programs in `c/` that link them (`run_detector`, `bench`, `throughput`) are therefore distributed under GPL-3.0 as built binaries. The MATLAB benchmark only calls ppg-beats and includes none of its code.
