# mzstab

Signal stability and time-window selection in time-resolved m/z spectra.

For each acquisition (about 30 s of scans, no chromatography), the program
follows the ions of the injected components over time, determines whether the
absolute intensity is stable, and recommends the time window that gives a
stable signal: globally across all acquisitions, and per target ion.

Thermo and Waters data are processed as separate, independent projects.

## Overview

Most users run the program from its graphical interface, which opens in the
web browser and runs entirely on the local computer. No file is edited by
hand: settings, components and samples are entered in forms and tables.

```
.raw files  ->  conversion to mzML  ->  mzstab interface  ->  report
```

## 1. Get the program

**Option A: standalone application (no Python needed)**

Download the zip for your system from the
[Releases](https://github.com/QDher/MSparser/releases) page, unzip it, and start
`mzstab-app.exe` (Windows) or `mzstab-app` (macOS, Linux). A console window
shows the log and the interface opens in the browser. Closing the console
stops the program.

**Option B: with Python (3.10 or later)**

```bash
pip install ".[ui]"
```

```bash
mzstab ui path/to/project
```

## 2. Install the converter

The program cannot read vendor `.raw` files directly. They are converted to
mzML by an external tool, which must be installed separately:

| Vendor | Converter | Notes |
|---|---|---|
| Thermo | [ThermoRawFileParser](https://github.com/compomics/ThermoRawFileParser) | Windows, Linux, macOS |
| Waters | msconvert ([ProteoWizard](https://proteowizard.sourceforge.io)) | Windows only for Waters files |

If mzML files already exist, the converter is not needed: place them in the
`mzml` folder of the project.

## 3. Use the interface

A *project folder* holds one instrument run. By default it is
`mzstab_project` in the home folder; another folder can be typed in the
sidebar. Put the `.raw` files in its `raw` subfolder (or mzML files in
`mzml`). Then follow the tabs from left to right:

| Tab | What to do |
|---|---|
| **1. Settings** | Enter the m/z windows (default ±10 and ±5 ppm). Everything else is read from the data files or uses a default; see below. |
| **2. Components** | Enter every injected molecule or ion: name, neutral formula, possible charge states (+1, +2, −1, −2) and charge mode. |
| **3. Samples** | Click *Add files found*, then fill in the reaction and the components injected in each acquisition (empty = all). |
| **4. Targets** | Generate the list of ions to follow and untick those not wanted. |
| **5. Run** | *Convert raw files*, then *Run analysis*. The log appears below the buttons. |
| **6. Results** | Recommended window, stable percentage, thresholds, traces per ion, flags, comparison of the m/z windows, report download. |

The sidebar shows which steps are complete.

### Components: charge mode

- **proton**: the molecule is charged by gain or loss of H⁺
  (`[M+H]+`, `[M+2H]2+`, `[M-H]-`, `[M-2H]2-`).
- **intrinsic**: the species is itself an ion (for example Na, or a quaternary
  ammonium). Its charge is applied by electron mass only, and it can form
  adducts with the other molecules (for example `[M+Na]+`).

Only entered species are used: no adduct is proposed with an ion that was not
entered. Dimers are not generated.

### Settings

Only the m/z windows have to be chosen. The other settings are automatic and
can be changed under *Advanced settings*:

| Setting | Automatic value |
|---|---|
| Vendor | Recognised from the files (Thermo or Waters) |
| Polarity, scan range | Read from the files |
| Window lengths | From the scan rate and the acquisition duration: lengths holding at least the minimum number of scans and no longer than 60 % of the shortest acquisition |
| Step between windows | 1 s, or the scan interval if longer |
| Minimum scans per window | 8 |
| Resolving power | Read from the file if present; otherwise 60,000 at m/z 200 (Orbitrap) or 20,000 (TOF). Only used to flag targets too close to be separated |
| Maximum injection time (Orbitrap) | The largest injection time observed in the file |
| Waters functions | Function 1 |

After the files have been read, the Settings tab lists what was found in each
file (vendor, polarity, scan range, number of scans, duration, scan rate,
resolving power and its origin). The values used by a run are saved in
`report/settings_used.json`.

Not available from the files, optional: the saturation intensity (TOF) and
the Waters lock-mass m/z (needed for conversion only).

### Recommended first run

Use one file, run the analysis, and compare one trace with FreeStyle or
MassLynx. Then add the other files. The validation step needs at least four
acquisitions.

## 4. Read the results

The **Results** tab shows, for each m/z window and time axis (acquisition
time, or time after signal onset):

- the **recommended global window** and the percentage of traces that are
  stable in it;
- the **thresholds derived from the data** (RSD limit and slope limit), with
  the figures to inspect them and fields to replace them manually;
- a **table per target ion**: its own window, compared with the global one;
- an **interactive chart of the traces** of a selected ion, with the stable
  window of each acquisition highlighted;
- **detection and flags** per ion: out of scan range, not detected, isotope
  mismatch, saturated, injection time at maximum, unresolved, multiple peaks;
- the **comparison of the m/z windows** (signal lost and interference removed
  by the narrow window);
- the **validation**: windows derived on half of the acquisitions and tested
  on the other half.

### Manual window

The automatic windows can be complemented by a window chosen by the user, in
the **Manual window** section of the Results tab: a start and end time (on
acquisition time) for all ions, or for single ions. A window set for one ion
takes precedence over the one set for all ions. *Apply and re-analyze*
evaluates the manual windows with the same criteria as the automatic ones.

The table below the form, also saved as `window_choice.csv`, lists per ion
the automatic window for all ions, the automatic window of that ion and the
manual window, each with the percentage of acquisitions in which the signal
is stable, the median RSD and slope, the mean intensity and its variation
between acquisitions. The column `chosen_by` marks each window as `automatic`
or `user`. Manual windows are drawn as dashed orange lines on the charts and
are recorded in `config.yaml`.

A manual window shorter than the automatic one can be classed as not stable
for weak ions: with few scans and a noisy signal, the absence of drift cannot
be established.

All tables and figures are also written to the `report` folder of the project
and can be downloaded as a zip from the same tab.

| File | Content |
|---|---|
| `trace_summary.csv` | Detection, mass error, isotope check and flags per trace |
| `trace_results.csv` | Phases, trend, selected window or `no_stable_window` per trace |
| `global_windows.csv`, `global_grid.csv` | Recommended window; stable percentage of every candidate window |
| `group_windows.csv` | Window per target ion, compared with the global window |
| `window_choice.csv` | Automatic and manual windows per ion, evaluated the same way; `chosen_by` marks the origin |
| `manual_window_traces.csv` | Per acquisition: mean, RSD, slope and stability in the manual windows |
| `window_comparison.csv` | Wide versus narrow m/z window |
| `thresholds.json` | Noise model, RSD factor, slope limit |
| `validation.csv`, `fixed_window.csv`, `mixed_model_*.txt` | Validation, reproducibility in the global window, dependence on intensity, m/z, reaction, run order |
| `*.png`, `traces_ppm*/` | Figures |

## How it works

1. **Targets** are computed from the entered components (monoisotopic m/z).
2. **Scans** are filtered: MS1 (and selected Waters functions), dominant
   polarity and scan range.
3. **Traces** are extracted per target: in each scan, the centroid closest to
   the theoretical m/z inside the window. Absolute intensities are used; there
   is no normalization. Undetected scans stay missing.
4. Each trace is **described**: rise, plateau and decay; mean, RSD and slope
   in sliding windows of several lengths.
5. **Thresholds** are derived from the data: the most stable window of each
   trace gives a noise model (`RSD² = a/I + b`) and the spread of slopes.
6. A window is **stable** if its RSD and slope are within the limits and its
   slope is estimated precisely enough. Per trace, the stable window with the
   lowest RSD/√n is selected. The global and per-ion windows are those stable
   in the highest percentage of traces.

## Files in a project folder

The interface writes plain-text files that record the inputs of a run. They
do not need to be edited by hand.

| File | Written by tab |
|---|---|
| `config.yaml` | Settings |
| `components.yaml` | Components |
| `samples.csv` | Samples |
| `work/targets.csv` | Targets |

## Advanced

### Command line

The same analysis runs without the interface, on the files above (examples in
`config/`):

```bash
mzstab convert -c config.yaml
```

```bash
mzstab targets -c config.yaml
```

```bash
mzstab run -c config.yaml
```

`run` executes `targets` (if no table exists), `ingest`, `extract` and
`analyze`; each can be called separately. `convert --dry-run` prints the
converter commands without running them. A project created in the interface
can be re-run with `mzstab run -c path/to/project/config.yaml`.

### CSV input

Spectra can be given as CSV instead of mzML: one file per acquisition in the
data folder, one row per centroid peak.

```
time,mz,intensity
0.20,195.0876,9876543
0.20,217.0696,987654
0.53,195.0877,9912345
```

Column names are detected (`time`/`rt`, `mz`/`m/z`, `intensity`, optional
`scan`). Time unit and polarity are set in the Settings tab. CSV files carry
no scan range or injection time, so the corresponding flags are unavailable.
An averaged spectrum exported from FreeStyle cannot be used: it has no time
dimension.

### Building the standalone application

An executable only runs on the system it was built on.

Windows:

```bash
packaging\build.bat
```

Linux or macOS:

```bash
bash packaging/build.sh
```

The result is the folder `dist/mzstab-app` (about 550 MB). The vendor
converters are not included.

Two workflows build the application on GitHub Actions, started from the
Actions tab:

| Workflow | Systems |
|---|---|
| `build-app` | Windows, Linux |
| `build-app-macos` | macOS, Apple silicon and Intel |

To attach the builds to a release, enter its version tag (for example `v1.0`)
in the *release_tag* field. macOS is separate because GitHub drops a macOS job
when no runner is available within 15 minutes; in that case, start
`build-app-macos` again later. On a Mac, `bash packaging/build.sh` builds the
application locally.

### Tests

```bash
pip install ".[test]"
```

```bash
python -m pytest -q
```

The tests generate synthetic mzML files (`tests/synth.py`) and run the full
pipeline. `python tests/synth.py demo` creates the same synthetic project,
which can be opened with `mzstab ui demo`.

## Status and limits

- Tested on synthetic data only; not yet validated on real Thermo or Waters
  files.
- The conversion step calls external tools and is not covered by the tests.
- Input must be centroided. Profile spectra produce a warning.
- Ion suppression and source instability are not corrected: the results
  describe when the signal is stable, not why it changes.
- Ions formed with species that were not entered (for example sodium from
  glassware) are not followed.
- The noise model is fitted on the lowest-RSD window of each trace, which
  biases it low; a multiplying factor compensates. Inspect the threshold
  figures before accepting the limits.
- Adducts combining charge carriers of opposite sign (for example
  `[M-2H+Na]-`) are not generated.
