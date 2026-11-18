# mzstab

Analysis of time-resolved m/z spectra acquired without chromatography
(direct infusion, ambient ionization, reaction monitoring; typically about
30 s of scans per acquisition), for Thermo and Waters instruments.

The program has two modes:

| Mode | Question answered |
|---|---|
| **Stability analysis** | In which time window of the acquisition is the signal stable? |
| **Intensity extraction** | What is the intensity of every possible ion in every file? |

In both modes you enter the molecules and ions that were injected; the
program computes the ions they can form (protonated forms, adducts, pairs)
and follows them in the data.

The program runs on your own computer and opens in the web browser. Nothing
is sent over the internet.

**Contents**

- [Windows: getting started](#windows-getting-started)
- [Using the program](#using-the-program)
- [Mode 1: stability analysis](#mode-1-stability-analysis)
- [Mode 2: intensity extraction](#mode-2-intensity-extraction)
- [macOS](#macos)
- [Linux](#linux)
- [Troubleshooting](#troubleshooting)
- [Reference](#reference)
- [For developers](#for-developers)
- [Status and limits](#status-and-limits)

---

## Windows: getting started

No Python and no command line are needed.

### 1. Download and start the program

1. Open the [Releases page](https://github.com/QDher/MSparser/releases) and
   choose the most recent version.
2. Under **Assets**, download `mzstab-app-windows-x64-<version>.zip`.
   Do not take *Source code (zip)*: it contains the program's source, not the
   program.
3. Right-click the downloaded zip and choose **Extract All**. The program
   does not work when started from inside the zip.
4. Open the extracted `mzstab-app` folder and double-click `mzstab-app.exe`
   (Windows may show it as `mzstab-app`, type *Application*). The `_internal`
   folder must stay next to it.
5. If Windows shows a blue *Windows protected your PC* message, click
   **More info**, then **Run anyway**. The program is not digitally signed.

A black console window opens, and after some seconds the program appears in
your web browser. The first start can take up to a minute. Keep the console
window open while working; closing it stops the program.

There is nothing to install or uninstall. To update, download the new zip
and use the new folder; the old one can be deleted. Your projects are kept,
because they are stored separately (see [Project folder](#project-folder)).

### 2. Install the file converter

The program cannot read the instruments' `.raw` files directly. They are
first converted to the open mzML format by a separate free tool. The program
runs the tool for you; it only has to be present on the computer.

**Thermo files: ThermoRawFileParser**

1. Open the
   [ThermoRawFileParser releases page](https://github.com/compomics/ThermoRawFileParser/releases).
2. Under **Assets** of the latest release, download the file ending in
   `-win.zip`.
3. Extract it to a permanent folder, for example
   `C:\Tools\ThermoRawFileParser\`.
4. If Windows blocks the program later (*Access denied*): right-click
   `ThermoRawFileParser.exe` > **Properties** > tick **Unblock** > **OK**.

**Waters files: msconvert**

Install [ProteoWizard](https://proteowizard.sourceforge.io), which contains
`msconvert.exe`.

**Tell the program where the converter is**

In the program, select the mode *Stability analysis*, open the tab
**1. Settings** > **Advanced settings** > **Converters**, and enter the full
path of the converter, for example
`C:\Tools\ThermoRawFileParser\ThermoRawFileParser.exe`. Giving the folder
that contains it also works. Click **Save settings**.

Tip: in File Explorer, hold Shift, right-click the `.exe`, and choose
**Copy as path** to get the exact path.

This setting is saved per project, so it is entered once for each new
project folder.

If your files are already in mzML format, no converter is needed.

---

## Using the program

### Project folder

The field at the top of the left panel is the **project folder**. It holds
the settings and the results of one study, not the data:

| In the project folder | Elsewhere, wherever you keep it |
|---|---|
| Settings, list of components, list of ions | The raw files |
| Intermediate files (`work`) | The converted mzML files |
| Result tables and figures (`report`) | |

By default it is `mzstab_project` in your user folder. To start another
study, type another folder name; it is created if needed. Your data is never
moved or copied: you indicate its folder inside the program.

### Choosing the mode

Below the project folder, **Mode** switches between *Stability analysis* and
*Intensity extraction*. Each mode has its own tabs, to be followed from left
to right. The list under the mode shows which steps are complete.

### Entering the components (both modes)

In the **Components** tab, add one row per molecule or ion that was injected:

| Column | What to enter |
|---|---|
| Name | A name of your choice |
| Neutral formula | The formula without charge: `C8H10N4O2`, `Na` (not `Na+`) |
| +1, +2, −1, −2 | Tick the charge states that are possible |
| The component is a… | *Molecule* or *Ion*, see below |

- **Molecule (gains or loses H⁺):** a neutral molecule that is protonated or
  deprotonated in the source, giving `[M+H]+`, `[M-H]-`, and so on.
- **Ion (already charged):** a species that is an ion by itself, such as Na,
  K, NH4, Cl, formate or a quaternary ammonium. No H⁺ is added or removed.
  **Choose this type for an ion to obtain its adducts with the molecules:**
  sodium with +1 gives `[M+Na]+`; chloride with −1 gives `[M+Cl]-`.

Examples:

| Name | Formula | Tick | Type |
|---|---|---|---|
| caffeine | `C8H10N4O2` | +1 | Molecule |
| sodium | `Na` | +1 | Ion |
| chloride | `Cl` | −1 | Ion |
| formate | `CHO2` | −1 | Ion |

Click **Save components**.

### The ions followed (both modes)

In the **Targets** tab, click **Generate targets**. The program proposes
three types of ions:

| Type | Content | Example |
|---|---|---|
| `ion` | One component, charged | `[M+H]+`, `Na+` |
| `adduct` | A molecule with one or two entered ions | `[M+Na]+`, `[M+H+Na]2+`, `[M+Cl]-` |
| `pair` | Two different components | `[M1+M2+H]+`, a cation with an anion such as `[Ca+Cl]+` |

Untick **use** for the ions you do not want, then click **Save target
selection**.

- Positive ions are followed in positive-mode files, negative ions in
  negative-mode files. The polarity is read from each file.
- Only the species you entered are used. No adduct is proposed with an ion
  that was not entered.
- A component is not combined with itself (no `[2M+H]+`), and ions made of
  three or more components are not proposed.
- After changing the components, generate the targets again.

---

## Mode 1: stability analysis

Finds the time window of the acquisition in which the signal is stable, for
all ions together and for each ion.

| Tab | What to do |
|---|---|
| **1. Settings** | Nothing has to be entered. Open *Advanced settings* only to give the converter path or to change a default. |
| **2. Components** | Enter the molecules and ions (see above). |
| **3. Samples** | Indicate the folders and list the files (see below). |
| **4. Targets** | Generate and review the ions (see above). |
| **5. Run** | Click **Convert raw files**, then **Run analysis**. |
| **6. Results** | Read the results (see below). |

### Samples tab: choosing the files

1. **Folder containing the raw files:** a subfolder of the project folder, or
   a full path such as `D:\data\experiment1`. The number of raw files found
   is shown below the field.
2. **Folder for the converted files:** where the mzML files are written. If
   you already have mzML files, give their folder here.
3. Click **Add files found in these folders**. Each file becomes a row.
4. Optionally fill in the other columns, then click **Save samples**.

Only the files directly in the folder are listed; subfolders are not read in
this mode. Rows can be deleted to leave files out.

| Column | Meaning |
|---|---|
| file | File name without extension |
| reaction | Optional label for the reaction or condition |
| components | Optional: names of the components injected in that file, separated by `;`. Empty means all components |
| run_order | Optional: position in the measurement sequence |

### Results tab

At the top, choose the m/z window (if several were set) and the time axis
(acquisition time, or time after the signal appears).

| Section | Content |
|---|---|
| **Global window** | The recommended time window and the percentage of traces that are stable in it. A *trace* is one ion in one file |
| **Scaled intensity profiles** | All traces, each divided by its own typical level. The ions shown can be chosen |
| **Traces with a stable signal in each window** | One cell per candidate window (start time × length); the colour is the percentage of traces stable in it |
| **Stability limits** | The two limits that define "stable", derived from the data and editable (see below) |
| **Manual window** | A window of your choice, evaluated next to the automatic ones |
| **Per target ion** | The window of each ion, and a chart of its traces |
| **Per acquisition** | For each file, how many ions are stable |
| **Detection and flags** | For each ion, in how many files it was detected or flagged |

**Stability limits.** A window is stable when its RSD is below the RSD limit
and its drift is within the slope limit. Both are derived from the data and
shown in two fields. To use your own values, change them and click **Apply
and re-analyze**; every table and exported file then uses the new limits. An
emptied field returns to the value derived from the data. On the *RSD versus
intensity* chart, hover or click a dot to see which file and ion it is.

**Manual window.** Enter a start and end time for all ions, or for single
ions. The table lists the automatic and the manual windows side by side, and
can be downloaded as `window_choice.csv`.

**Download.** *Prepare report archive* at the bottom gives all tables and
figures as one zip. They are also in the `report` folder of the project.

Recommended first run: use one file and compare one trace with FreeStyle or
MassLynx before adding the others.

---

## Mode 2: intensity extraction

Reports the intensity of every possible ion in every file, whatever its
stability. Unstable signals are flagged, not removed. This mode is meant for
many files (up to about a hundred) organised in folders.

| Tab | What to do |
|---|---|
| **1. Components** | Enter the molecules and ions (see above). |
| **2. Files** | Indicate the main folders and tick the folders to use (see below). |
| **3. Targets** | Generate and review the ions (see above). |
| **4. Extract** | Click **Extract intensities**. |
| **5. Results** | Read and download the intensities (see below). |

### Files tab: choosing the folders

1. **Folder containing the raw files:** the top folder of your data. Its
   subfolders, at any depth, are found automatically.
2. **Folder for the converted files:** where the mzML files are written. The
   folder tree of the raw files is reproduced there.
3. In the table **Folders to use**, tick the folders wanted. Each row shows
   the number of raw and converted files in that folder.
4. Click **Convert the raw files of the ticked folders**.

- Nothing is used until it is ticked. *Tick all* selects everything.
- A ticked folder does not include its subfolders; tick each one wanted.
- Every file of a ticked folder is used. The same components apply to all
  files.

### How an intensity is obtained

For each ion in each file:

1. **Adaptive m/z window.** The windows ±5, 10, 15 and 20 ppm are tried in
   this order. The first one containing a signal is used and reported. A
   signal is a peak present in at least 5 consecutive scans (1, 3 or 5,
   selectable).
2. **Peaks in the window.** If the window holds several peaks, the one
   closest to the theoretical m/z, the most intense one, and any peak of at
   least 10 % of the most intense are listed.
3. **Intensity.** The acquisition is covered by six time windows, each
   lasting one third of the acquisition and spread evenly from the first to
   the last scan (for 30 s: six 10 s windows starting every 4 s). The
   intensity is the mean over the window with the lowest RSD.
4. **Flags.** Each value comes with its RSD, the spread between the six
   windows, the m/z window used and the mass error.

| Flag | Meaning |
|---|---|
| `low_stability` | RSD in the window above 30 % |
| `window_variation` | The six windows differ by 25 % or more of the reported intensity |
| `gaps` | Signal missing in more than 10 % of the scans of the window |
| `multiple_peaks` | Several peaks in the m/z window |
| `wide_mz_window` | Found only with a window wider than the first one |
| `isotope_mismatch` | Isotope peak ratio differs from the theoretical one |
| `saturated` | Intensity at the saturation level (if one was set) |
| `short_acquisition` | Too few scans per window; all scans of the acquisition were used |

The numbers in steps 1, 2 and 4 can be changed in the **Extract** tab.

### Results tab

| Section | Content |
|---|---|
| **Intensity table** | *files × ions*: one row per file, one column per ion. *One row per value*: with flags, RSD, variation, m/z window, mass error. Filters by folder and ion; download button |
| **Choice of the peak** | When several peaks were found: a general rule (closest, or most intense) and a table to tick the peak to export for single ions and files |
| **Detail of one ion in one file** | The signal over the acquisition with the window used, and the mean of each of the six windows |

Files already processed with the same settings are reused, so adding files
and extracting again only processes the new ones. A file that cannot be read
is listed with its error and does not stop the others.

---

## macOS

There is currently no ready-made macOS application on the Releases page.
Use one of the two options below; both need the Terminal.

**Option A: run with Python** (3.10 or later)

1. Download *Source code (zip)* from the
   [Releases page](https://github.com/QDher/MSparser/releases) and unzip it.
2. In the Terminal, go to the unzipped folder and install:

   ```bash
   pip install ".[ui]"
   ```

3. Start the program on a project folder of your choice:

   ```bash
   mzstab ui ~/mzstab_project
   ```

**Option B: build the application on your Mac**

In the unzipped source folder:

```bash
bash packaging/build.sh
```

The application is created in `dist/mzstab-app`; start `mzstab-app` inside
it. macOS may ask for confirmation the first time (right-click > **Open**).

**Converters on macOS**

- **Thermo:** download the ThermoRawFileParser file ending in `-osx.zip`
  (Intel) or `-osx-arm64.zip` (Apple silicon) from its
  [releases page](https://github.com/compomics/ThermoRawFileParser/releases),
  unzip it, and give the path of the `ThermoRawFileParser` program under
  *Settings > Advanced settings > Converters*.
- **Waters:** conversion is not possible on macOS. Convert the files to mzML
  on a Windows computer and give the folder of the mzML files in the program.

Using the program is then identical to Windows.

---

## Linux

**Application without Python**

1. Download `mzstab-app-linux-x64-<version>.zip` from the
   [Releases page](https://github.com/QDher/MSparser/releases) and unzip it.
2. Start it from a terminal:

   ```bash
   chmod +x mzstab-app/mzstab-app
   ```

   ```bash
   ./mzstab-app/mzstab-app
   ```

The program opens in the browser; if it does not, open the address shown in
the terminal (usually `http://localhost:8501`). A project folder can be given
as argument: `./mzstab-app/mzstab-app ~/my_project`.

**With Python** (3.10 or later), from the unzipped source code:

```bash
pip install ".[ui]"
```

```bash
mzstab ui ~/mzstab_project
```

**Converters on Linux**

- **Thermo:** download the ThermoRawFileParser file ending in `-linux.zip`
  from its
  [releases page](https://github.com/compomics/ThermoRawFileParser/releases),
  unzip it, and give the path of the `ThermoRawFileParser` program under
  *Settings > Advanced settings > Converters*.
- **Waters:** conversion is not possible natively on Linux. Convert the files
  to mzML on a Windows computer.

Using the program is then identical to Windows.

---

## Troubleshooting

| Problem | Cause and solution |
|---|---|
| The unzipped folder contains `README.md` and `pyproject.toml` but no `.exe` | The *Source code* zip was downloaded. Take `mzstab-app-windows-x64-<version>.zip` |
| Nothing happens when starting the program | Extract the whole zip first. Wait up to a minute at the first start. Check that the antivirus does not block it |
| *Windows protected your PC* | Click **More info** > **Run anyway** |
| The browser page does not open | Open `http://localhost:8501` in the browser while the console window is open |
| *Convert failed: no .raw file in …* | The raw folder is empty or wrong. Check the folder field; the path used is shown in the Run tab |
| *converter not found* | Enter the full path of the converter under *Settings > Advanced settings > Converters* |
| *Access denied* (*Accès refusé*) at conversion | The converter setting is not the program file, or Windows blocks it: right-click the `.exe` > **Properties** > **Unblock** |
| No adduct such as `[M+Na]+` in the targets | The ion must be entered with the type *Ion (already charged)*; then generate the targets again |
| An ion is "not detected" | Check its formula and charge, the polarity of the file and the scan range. In the stability mode, see *Detection and flags*; a mass error above the m/z window also gives this result |
| An ion is "out of scan range" | Its m/z is outside the range acquired, for example Na⁺ at m/z 23. Its adducts are still followed |
| "No stable window" for every trace | The stability limits are too strict for the data, or the signal drifts. Adjust the limits in the Results tab |
| Results look outdated after changing a setting | Run the analysis or extraction again; in the Targets tab, regenerate the targets after changing components |

If a step fails, the message shown in the program states the reason. Include
it, or a screenshot, when reporting a problem.

---

## Reference

### Automatic settings (stability mode)

No setting has to be entered. These values are used unless changed under
*Advanced settings*; an empty field means automatic.

| Setting | Automatic value |
|---|---|
| m/z window | ±15 ppm around the theoretical m/z. Several values (for example `10, 5`) are all run and compared |
| Vendor, polarity, scan range | Read from the files |
| Time-window lengths | From the scan rate and the acquisition duration: lengths holding at least the minimum number of scans and no longer than 60 % of the shortest acquisition |
| Step between windows | 1 s, or the scan interval if longer |
| Minimum scans per window | 8 |
| Resolving power | Read from the file if present; otherwise 60,000 at m/z 200 (Orbitrap) or 20,000 (TOF). Only used to flag ions too close to be separated |
| Maximum injection time (Orbitrap) | The largest value observed in the file |
| Waters functions | Function 1 |

After the files are read, the Settings tab lists what was found in each file.
Optional and not available from the files: the saturation intensity (TOF)
and the Waters lock-mass m/z (used at conversion).

### How the stability analysis works

1. The m/z of each ion is computed from the entered components.
2. Scans are filtered: MS1 (and the selected Waters functions), dominant
   polarity and scan range.
3. For each ion and file, a trace is extracted: in each scan, the peak
   closest to the theoretical m/z inside the window. Absolute intensities are
   used, without normalization.
4. Each trace is described: rise, plateau and decay; mean, RSD and slope in
   sliding windows of several lengths.
5. The most stable window of each trace gives a noise model (the typical RSD
   at a given intensity). From it, the RSD limit and the slope limit are
   derived, unless set by the user.
6. A window is stable if its RSD and slope are within the limits and the
   slope is measured precisely enough. Per trace, the stable window with the
   lowest uncertainty of the mean is selected. The global and per-ion windows
   are those stable in the highest percentage of traces.

### Output files

All are in the `report` folder of the project.

Stability analysis:

| File | Content |
|---|---|
| `global_windows.csv`, `global_grid.csv` | Recommended window; stable percentage of every candidate window |
| `group_windows.csv` | Window per ion, compared with the global window |
| `window_choice.csv` | Automatic and manual windows per ion; `chosen_by` marks the origin |
| `trace_results.csv` | Per ion and file: phases, trend, selected window or `no_stable_window` |
| `file_summary.csv` | Per file: analyzed ions, ions with a stable window |
| `trace_summary.csv` | Per ion and file: detection, mass error, isotope check, flags |
| `thresholds.json`, `threshold_points.csv` | Limits used and the data behind the RSD chart |
| `manual_window_traces.csv` | Per ion and file, values in the manual windows |
| `window_comparison.csv` | Comparison of the m/z windows, when several are set |
| `validation.csv`, `fixed_window.csv`, `mixed_model_*.txt` | Validation on held-out files, reproducibility, dependence on intensity, m/z, reaction and run order |
| `settings_used.json` | The settings used by the run |
| `*.png`, `traces_ppm*/` | Figures |

Intensity extraction:

| File | Content |
|---|---|
| `intensity_matrix.csv` | Files as rows, ions as columns, intensities as values |
| `intensity_values.csv` | One row per file and ion: intensity, flags, RSD, variation, m/z window, mass error, window used; undetected ions are listed with their status |
| `intensity_peaks.csv` | All candidate peaks, with the one selected |
| `intensity_window_means.csv` | Mean, RSD and slope of each of the six windows of every reported peak |
| `intensity_files.csv` | Per file: folder, processing status or error, polarity, scan range, scan rate |

### Files written in the project folder

The program records the inputs of a run in plain-text files. They do not
need to be edited by hand.

| File | Written by |
|---|---|
| `config.yaml` | Settings, folders, limits, mode |
| `components.yaml` | Components tab |
| `samples.csv` | Samples tab |
| `work/targets.csv` | Targets tab |

### CSV input

Spectra can be given as CSV instead of mzML: one file per acquisition, one
row per centroid peak.

```
time,mz,intensity
0.20,195.0876,9876543
0.20,217.0696,987654
0.53,195.0877,9912345
```

Column names are detected (`time`/`rt`, `mz`/`m/z`, `intensity`, optional
`scan`). The time unit and the polarity are set under *Advanced settings*.
CSV files carry no scan range or injection time, so the corresponding flags
are unavailable. An averaged spectrum exported from FreeStyle cannot be
used: it has no time dimension.

---

## For developers

### Install from source

```bash
pip install ".[ui,test]"
```

### Command line

Both modes run without the interface, on the files of a project folder
(examples in `config/`):

```bash
mzstab convert -c config.yaml
```

```bash
mzstab targets -c config.yaml
```

```bash
mzstab run -c config.yaml
```

```bash
mzstab intensity -c config.yaml --jobs 4
```

`run` executes the stability analysis (`targets` if no table exists, then
`ingest`, `extract`, `analyze`; each can be called separately). `intensity`
runs the intensity extraction and needs `mode: intensity` in the
configuration for recursive conversion; `--jobs` processes files in
parallel. `convert --dry-run` prints the converter commands without running
them. `mzstab ui <project>` starts the interface.

In `components.yaml`, the two component types are stored as `mode: proton`
(molecule) and `mode: intrinsic` (ion); `molecule` and `ion` are accepted.

### Tests

```bash
python -m pytest -q
```

The tests generate synthetic mzML files (`tests/synth.py`) and run both
modes, including the interface. `python tests/synth.py demo` creates a
synthetic project that can be opened with `mzstab ui demo`.

### Building the application

An executable only runs on the system it was built on.

Windows:

```bash
packaging\build.bat
```

Linux or macOS:

```bash
bash packaging/build.sh
```

The result is the folder `dist/mzstab-app` (about 550 MB). The converters are
not included.

Two GitHub Actions workflows build the application, started from the Actions
tab:

| Workflow | Systems |
|---|---|
| `build-app` | Windows, Linux |
| `build-app-macos` | macOS, Apple silicon and Intel |

To publish a release, enter a version tag (for example `v1.4`) in the
*release_tag* field: the builds are attached to the release of that tag,
which is created as a pre-release if needed. macOS is separate because GitHub
drops a macOS job when no runner is available within 15 minutes; so far no
macOS build has obtained a runner.

---

## Status and limits

- Developed and tested on synthetic data. Use on real data is recent and
  limited to a few Thermo files; Waters files and the macOS application have
  not been tried.
- The conversion step calls external tools and is not covered by the tests.
- Input must be centroided. Profile spectra produce a warning.
- Ion suppression and source instability are not corrected: the results
  describe the signal, not why it changes.
- Ions formed with species that were not entered (for example sodium from
  glassware) are not followed.
- Ions combining charge carriers of opposite sign (for example
  `[M-2H+Na]-`), a component with itself, or three or more components are not
  generated.
- Stability mode: with few traces, the limits derived from the data are not
  reliable; the program warns and the limits can be set by hand.
- Intensity mode: the same components apply to all files; files are
  processed one after another in the application.
