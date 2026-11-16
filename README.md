# mzstab

Signal stability and time-window selection in time-resolved m/z spectra
(targeted, absolute intensities, one independent run per vendor).

## Install

```bash
pip install -e ".[test]"
```

## Inputs (see `config/`)

| File | Content |
|---|---|
| `thermo.yaml` / `waters.yaml` | One configuration per vendor run |
| `components.yaml` | Injected components: neutral formula, allowed charges, charge mode |
| `samples.csv` | `file` (name without extension of the mzML or CSV), `reaction`, `components` (`;`-separated, empty = all), `run_order` |

## Run

```bash
mzstab convert -c config/thermo.yaml
```

```bash
mzstab targets -c config/thermo.yaml
```

Review `work/targets.csv` and set `use` to 0 for ions not to follow. Then:

```bash
mzstab run -c config/thermo.yaml
```

`run` executes `targets` (if no table exists), `ingest`, `extract` and `analyze`.
Each step can also be called separately. `convert --dry-run` prints the
converter commands without running them.

## Graphical interface

```bash
pip install -e ".[ui]"
```

```bash
mzstab ui path/to/project
```

This opens a page in the browser, running locally. A project folder holds one
vendor run (`config.yaml`, `components.yaml`, `samples.csv`, `mzml/`, `work/`,
`report/`); it is created if needed. The tabs follow the workflow:
Settings, Components, Samples, Targets, Run, Results. The interface writes the
same files as the command line, so a run can be repeated with
`mzstab run -c path/to/project/config.yaml`.

In the Results tab, the derived thresholds can be replaced by manual values
(stored under `thresholds: override:` in `config.yaml`; the derived values are
kept in `thresholds.json`).

## Standalone application

The interface can be packaged as a program that needs no Python installation.
An executable only runs on the operating system it was built on, so build it
on the target system.

Windows (from the repository root):

```bash
packaging\build.bat
```

Linux or macOS:

```bash
bash packaging/build.sh
```

The result is the folder `dist/mzstab-app` (about 550 MB; about 200 MB zipped).
Distribute the whole folder. Start `mzstab-app` (`mzstab-app.exe` on Windows);
a console window shows the log and the page opens in the browser. Closing the
console stops the program. The project folder defaults to
`<home>/mzstab_project`; another folder can be given as argument or typed in
the sidebar.

`.github/workflows/build-app.yml` builds the three systems on GitHub Actions.
The vendor converters are not included in the application.

## CSV input

Spectra can be given as CSV instead of mzML: one file per acquisition in the
data folder, one row per centroid peak.

```
time,mz,intensity
0.20,195.0876,9876543
0.20,217.0696,987654
0.53,195.0877,9912345
```

Column names are detected (`time`/`rt`, `mz`/`m/z`, `intensity`, optional
`scan`); set `csv: time_col:` etc. in the configuration if detection fails.
`csv: time_unit:` is `s` or `min`; `csv: polarity:` is `1` or `-1`. Rows are
grouped into scans by the scan column, or by identical time values. CSV files
carry no scan range or injection time, so the out-of-range and injection-time
flags are not available, and all rows are treated as MS1.

## Method

1. **Targets.** `[M+H]+`, `[M+2H]2+`, `[M-H]-`, `[M-2H]2-` for the allowed
   charge states; intrinsic ions by electron mass only; adducts only with
   entered intrinsic ions (same sign as the total charge). No dimers.
2. **Scan filter.** Configured MS level and Waters functions, then the dominant
   polarity and scan range. Dropped scans are counted in `scans_dropped.csv`.
3. **Extraction.** Per scan, the centroid closest to the theoretical m/z inside
   each ppm window (with an absolute floor). Missing values stay missing.
4. **Flags.** Not detected, out of scan range, saturated, injection time at
   maximum, unresolved pair of targets, multiple centroids, isotope mismatch.
   A trace is analyzed if detected in >= `min_consecutive_scans` consecutive
   scans and its isotope ratio is not contradicted.
5. **Description.** Rise/plateau/decay by change-point detection; sliding
   windows (mean, RSD, relative slope and its standard error) on acquisition
   time and on time after signal onset; linear trend over the acquisition.
6. **Thresholds (from the data).** At the reference window length, the most
   stable window of each trace is used to fit `RSD^2 = a/I + b` and to measure
   the spread of slopes. A window is stable if its RSD is below
   `rsd_factor x model`, its |slope| is below `slope_limit`, and the slope
   standard error is also below `slope_limit`.
7. **Selection.** Per trace: the stable window with the lowest RSD/sqrt(n).
   Global and per target ion: the window that is stable in the highest
   percentage of traces.
8. **Validation.** Thresholds and windows derived on a random half of the
   acquisitions are evaluated on the other half.

## Outputs (`report_dir`)

| File | Content |
|---|---|
| `trace_summary.csv` | Detection, mass error, isotope check and flags per trace and ppm setting |
| `trace_results.csv` | Phases, trend, selected window or `no_stable_window`, per trace and time axis |
| `thresholds.json` | Noise model, RSD factor, slope limit |
| `global_windows.csv`, `global_grid.csv` | Recommended global window; stable percentage for every candidate window |
| `group_windows.csv` | Window per target ion and comparison with the global window |
| `window_comparison.csv` | Wide versus narrow m/z window: signal lost, interference removed, effect on stability |
| `validation.csv` | Train/test percentages |
| `fixed_window.csv` | Between-acquisition CV of the mean intensity in the global window |
| `mixed_model_*.txt` | Dependence on intensity, m/z, reaction, run order |
| `*.png`, `traces_ppm*/` | Profiles, stable-percentage heatmaps, RSD versus intensity, threshold histograms, mass error, traces per target |

Results are produced for every ppm setting (`ppm` column or `_ppmN` suffix).

## Tests

```bash
python -m pytest -q
```

The tests generate synthetic centroided mzML files (`tests/synth.py`) and run
the full pipeline. To create the same synthetic project for inspection:

```bash
python tests/synth.py demo && mzstab run -c demo/config.yaml
```

## Status and limits

- Tested on synthetic mzML only. Not yet run on real Thermo or Waters files.
- The converter commands (`convert`) are not tested here; check them with
  `--dry-run` against the installed ThermoRawFileParser / msconvert versions.
- Input must be centroided. Profile spectra produce a warning.
- The noise model is fitted on the lowest-RSD window of each trace, which
  biases it low; `rsd_factor` compensates. Inspect
  `threshold_histograms_*.png` and `rsd_vs_intensity_*.png` before accepting
  the thresholds.
- Without `max_injection_time_ms`, the injection-time flag uses the largest
  observed value as the maximum and is only indicative.
- Adducts combining carriers of opposite sign (for example `[M-2H+Na]-`) are
  not generated.
