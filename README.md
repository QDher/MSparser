# MSparser

Analysis of time-resolved m/z spectra acquired without chromatography
(direct infusion, ambient ionization, reaction monitoring; typically about
30 s of scans per acquisition), for Thermo and Waters instruments.

The program has two modes:

| Mode | Question answered |
|---|---|
| **Intensity extraction** (default) | What is the intensity of every possible ion in every file? |
| **Stability analysis** | In which time window of the acquisition is the signal stable? |

In both modes you enter the molecules and ions that were injected; the
program computes the ions they can form (protonated forms, adducts, pairs)
and follows them in the data.

The program runs on your own computer and opens in the web browser. Nothing
is sent over the internet.

**Contents**

- [Windows: getting started](#windows-getting-started)
- [Using the program](#using-the-program)
- [Intensity extraction](#intensity-extraction)
- [Stability analysis](#stability-analysis)
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
2. Under **Assets**, download `MSparser-windows-x64-<version>.zip`.
   (Versions up to v1.3 were named *mzstab*: their file is
   `mzstab-app-windows-x64-<version>.zip` and the program `mzstab-app.exe`.)
   Do not take *Source code (zip)*: it contains the program's source, not the
   program.
3. Right-click the downloaded zip and choose **Extract All**. The program
   does not work when started from inside the zip.
4. Open the extracted `MSparser` folder and double-click `MSparser.exe`
   (Windows may show it as `MSparser`, type *Application*). The `_internal`
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

**Let the program find the converter**

The simplest way: put the converter's folder inside the folder of the
program, next to `MSparser.exe`. It is then found automatically, in every
project.

```
MSparser\
  MSparser.exe
  _internal\
  ThermoRawFileParser\        <- the unzipped converter
    ThermoRawFileParser.exe
```

Otherwise, give its location in the program: open the **Converter** panel
(intensity mode: tab **1. Folders**; stability mode: tab **5. Run**) and enter
the program file or the folder that contains it, for example
`C:\Tools\ThermoRawFileParser`. This is saved with the settings of the data
folder or project.

The panel shows whether the converter was found and where. Tip: in File
Explorer, hold Shift, right-click the folder or the `.exe`, and choose
**Copy as path**.

When updating MSparser to a new version, move the converter folder into the
new `MSparser` folder, or keep it elsewhere and give its location.

If your files are already in mzML format, no converter is needed.

---

## Using the program

### Choosing the mode

At the top of the left panel, **Mode** switches between *Intensity
extraction* (selected at first start) and *Stability analysis*. Each mode has
its own tabs, to be followed from left to right.

### Where your data is, and where the program writes

**Intensity extraction** works directly in your data folders. In the left
panel, give the **Data folder** that holds your acquisitions and tick
**Include subfolders** to process its subfolders as well. In every folder
that holds acquisitions, the program creates one subfolder named `MSparser`
and writes everything there:

```
C:\example_data\root\
  sample1.RAW
  sample2.RAW
  MSparser\                    <- created by the program
    mzml\                      converted files
    components.yaml            compound list of this folder
    config.yaml                settings used
    work\                      ion list and intermediate files
    results\                   intensity tables of this folder
    results_all\               tables combining root and its subfolders
  subfolder\
    sample3.RAW
    MSparser\                  the same, for this subfolder
```

Your raw files are not moved or modified. Deleting an `MSparser` folder
removes everything the program made for that folder. Recently used folders
are offered in the left panel.

**Stability analysis** uses a **project folder** instead, which holds the
settings and results of one study; the folders of the raw and converted
files are indicated inside the program. By default it is `MSparser_project`
in your user folder.

### Entering the compounds (both modes)

Add one row per molecule or ion that was injected:

| Column | What to enter |
|---|---|
| Name | A name of your choice |
| Neutral formula | The formula without charge: `C8H10N4O2`, `Na` (not `Na+`) |
| Charges | A list to pick from: every charge state the component can take, from −3 to +3 |
| The component is a… | *Molecule* or *Ion*, see below |

- **Molecule (gains or loses H⁺):** a neutral molecule that is protonated or
  deprotonated in the source, giving `[M+H]+`, `[M-H]-`, and so on.
- **Ion (already charged):** a species that is an ion by itself, such as Na,
  K, NH4, Cl, formate or a quaternary ammonium. No H⁺ is added or removed.
  **Choose this type for an ion to obtain its adducts with the molecules:**
  sodium with +1 gives `[M+Na]+` and `[M-2H+Na]-`; chloride with −1 gives
  `[M+Cl]-`.

Examples:

| Name | Formula | Charges | Type |
|---|---|---|---|
| caffeine | `C8H10N4O2` | +1 | Molecule |
| citric acid | `C6H8O7` | −1, −2, −3 | Molecule |
| sodium | `Na` | +1 | Ion |
| chloride | `Cl` | −1 | Ion |

Click **Save components**.

### The ions followed (both modes)

From the compounds, the program proposes the ions to follow:

| Type | Content | Examples |
|---|---|---|
| `ion` | One compound, charged by protons, or an entered ion itself | `[M+H]+`, `[M-2H]2-`, `Na+` |
| `adduct` | One molecule with protons gained or lost and entered ions attached | `[M+Na]+`, `[M+H+Na]2+`, `[M-2H+Na]-`, `[M-H+2Na]+`, `[M+Cl]-` |
| `pair` | Two different molecules, with protons and entered ions; or a cation with an anion | `[M1+M2+H]+`, `[M1+M2-2H+Na]-`, `[Ca+Cl]+` |

Rules:

- The total charge of an ion must be one of the charges picked for the
  molecule (for a pair: a charge of either molecule, or the sum of one charge
  of each).
- An ion holds at most two molecules, and never the same molecule twice (no
  `[2M+H]+`). Protons gained or lost and attached ions do not count as
  molecules.
- Up to 3 protons are gained or lost, and up to 3 entered ions attached.
  Cations and anions are not attached to the same ion.
- Only the species you entered are used. No adduct is proposed with an ion
  that was not entered.
- Positive ions are followed in positive-mode files, negative ions in
  negative-mode files. The polarity is read from each file.

The list can be long. In the ions tab, untick **use** for the ions you do not
want, then save the selection. After changing the compounds, the ions are
generated again.

---

## Intensity extraction

Reports the intensity of every possible ion in every file, whatever its
stability. Unstable signals are flagged, not removed. This mode is meant for
many files (up to about a hundred) organised in folders.

| Tab | What to do |
|---|---|
| **1. Folders** | Check the folders found, untick those to leave out, and check the converter. |
| **2. Compounds** | Enter the compound list of each folder. |
| **3. Ions** | Optional: review the ions of each folder and untick some. |
| **4. Convert and extract** | Click **Convert and extract**. |
| **5. Results** | Read and download the intensities, for all folders or one. |

### Folders

After the data folder is given in the left panel, tab 1 lists every folder
that holds acquisitions (raw, mzML or CSV files), with its numbers of files
and compounds. Untick **use** to leave a folder out. Each folder is processed
on its own and keeps its own files (see above).

### Compound lists

Each folder has its own compound list. In tab 2, choose a folder, enter its
compounds and save. **Copy this list to the other folders** gives the same
list to all of them; lists can then be changed folder by folder. A folder
without a compound list is skipped at extraction.

### Convert and extract

**Convert and extract** converts the raw files not yet converted, generates
the ions of each folder, and extracts the intensities. **Extract only** skips
the conversion. A table then gives, per folder, the outcome and the numbers
of files and values. Files already processed with the same settings and ions
are reused, so running again after adding files only processes the new ones.
The settings of tab 4 apply to all folders.

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

Choose **All folders** or one folder.

| Scope | Content |
|---|---|
| **All folders** | The combined intensity table (*files × ions*, or *one row per value* with flags), with filters by folder and ion, and a download button |
| **One folder** | Its intensity table; the **choice of the peak** when several peaks were found (a general rule, closest or most intense, and a table to tick the peak to export for single ions and files); the **detail of one ion in one file** (signal over the acquisition with the window used, and the mean of each of the six windows) |

A file that cannot be read is listed with its error and does not stop the
others.

---

## Stability analysis

Finds the time window of the acquisition in which the signal is stable, for
all ions together and for each ion.

| Tab | What to do |
|---|---|
| **1. Settings** | Nothing has to be entered. Open *Advanced settings* only to change a default. |
| **2. Components** | Enter the molecules and ions (see above). |
| **3. Samples** | Indicate the folders and list the files (see below). |
| **4. Targets** | Generate and review the ions (see above). |
| **5. Run** | Check the **Converter** panel, click **Convert raw files**, then **Run analysis**. |
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
   msparser ui ~/MSparser_project
   ```

**Option B: build the application on your Mac**

In the unzipped source folder:

```bash
bash packaging/build.sh
```

The application is created in `dist/MSparser`; start `MSparser` inside
it. macOS may ask for confirmation the first time (right-click > **Open**).

**Converters on macOS**

- **Thermo:** download the ThermoRawFileParser file ending in `-osx.zip`
  (Intel) or `-osx-arm64.zip` (Apple silicon) from its
  [releases page](https://github.com/compomics/ThermoRawFileParser/releases),
  unzip it, and give its folder in the **Converter** panel of the program.
- **Waters:** conversion is not possible on macOS. Convert the files to mzML
  on a Windows computer and give the folder of the mzML files in the program.

Using the program is then identical to Windows.

---

## Linux

**Application without Python**

1. Download `MSparser-linux-x64-<version>.zip` from the
   [Releases page](https://github.com/QDher/MSparser/releases) and unzip it.
2. Start it from a terminal:

   ```bash
   chmod +x MSparser/MSparser
   ```

   ```bash
   ./MSparser/MSparser
   ```

The program opens in the browser; if it does not, open the address shown in
the terminal (usually `http://localhost:8501`). A project folder can be given
as argument: `./MSparser/MSparser ~/my_project`.

**With Python** (3.10 or later), from the unzipped source code:

```bash
pip install ".[ui]"
```

```bash
msparser ui ~/MSparser_project
```

**Converters on Linux**

- **Thermo:** download the ThermoRawFileParser file ending in `-linux.zip`
  from its
  [releases page](https://github.com/compomics/ThermoRawFileParser/releases),
  unzip it, and give its folder in the **Converter** panel of the program.
- **Waters:** conversion is not possible natively on Linux. Convert the files
  to mzML on a Windows computer.

Using the program is then identical to Windows.

---

## Troubleshooting

| Problem | Cause and solution |
|---|---|
| The unzipped folder contains `README.md` and `pyproject.toml` but no `.exe` | The *Source code* zip was downloaded. Take `MSparser-windows-x64-<version>.zip` |
| Nothing happens when starting the program | Extract the whole zip first. Wait up to a minute at the first start. Check that the antivirus does not block it |
| *Windows protected your PC* | Click **More info** > **Run anyway** |
| The browser page does not open | Open `http://localhost:8501` in the browser while the console window is open |
| *Convert failed: no .raw file in …* | The raw folder is empty or wrong. Check the folder field; the path used is shown in the Run tab |
| *converter not found* | Put the converter's folder inside the `MSparser` folder, or give its location in the **Converter** panel |
| *Access denied* (*Accès refusé*) at conversion | The converter setting is not the program file, or Windows blocks it: right-click the `.exe` > **Properties** > **Unblock** |
| No adduct such as `[M+Na]+` in the ions | The ion must be entered with the type *Ion (already charged)*, and the total charge of the adduct must be among the charges picked for the molecule |
| A folder is "skipped" at extraction | It has no compound list: enter one in tab 2, or copy the list of another folder |
| *Permission denied* when extracting | The program cannot write in the data folder (read-only location). Copy the data to a folder where you can write |
| An ion is "not detected" | Check its formula and charge, the polarity of the file and the scan range. In the stability mode, see *Detection and flags*; a mass error above the m/z window also gives this result |
| An ion is "out of scan range" | Its m/z is outside the range acquired, for example Na⁺ at m/z 23. Its adducts are still followed |
| "No stable window" for every trace | The stability limits are too strict for the data, or the signal drifts. Adjust the limits in the Results tab |
| Results look outdated after changing a setting | Run the analysis or extraction again. In the stability mode, regenerate the targets after changing components |

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

Intensity extraction, in `<data folder>\MSparser\results` of each folder:

| File | Content |
|---|---|
| `intensity_matrix.csv` | Files as rows, ions as columns, intensities as values |
| `intensity_values.csv` | One row per file and ion: intensity, flags, RSD, variation, m/z window, mass error, window used; undetected ions are listed with their status |
| `intensity_peaks.csv` | All candidate peaks, with the one selected |
| `intensity_window_means.csv` | Mean, RSD and slope of each of the six windows of every reported peak |
| `intensity_files.csv` | Per file: processing status or error, polarity, scan range, scan rate |

And in `<top data folder>\MSparser\results_all`, for all folders together:
`intensity_matrix_all.csv`, `intensity_values_all.csv` and
`intensity_files_all.csv`, with a `data_folder` column.

### Files written in the project folder (stability analysis)

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
msparser convert -c config.yaml
```

```bash
msparser targets -c config.yaml
```

```bash
msparser run -c config.yaml
```

```bash
msparser folder C:\example_data\root --components compounds.yaml
```

```bash
msparser intensity -c config.yaml --jobs 4
```

- `run` executes the stability analysis on a project folder (`targets` if no
  table exists, then `ingest`, `extract`, `analyze`; each can be called
  separately).
- `folder` runs the intensity extraction on a data folder and its subfolders,
  as the interface does: each folder gets its `MSparser` subfolder.
  `--components` copies one compound list to every folder; without it, each
  folder uses the list already in its `MSparser` subfolder. `--no-subfolders`
  and `--no-convert` are available.
- `intensity` runs the intensity extraction on one configuration file, for
  example `<folder>/MSparser/config.yaml`, to process a single folder;
  `--jobs` processes files in parallel.
- `convert --dry-run` prints the converter commands without running them.
- `msparser ui [folder]` starts the interface.

In `components.yaml`, the two component types are stored as `mode: proton`
(molecule) and `mode: intrinsic` (ion); `molecule` and `ion` are accepted.

### Tests

```bash
python -m pytest -q
```

The tests generate synthetic mzML files (`tests/synth.py`) and run both
modes, including the interface. `python tests/synth.py demo` creates a
synthetic project that can be opened with `msparser ui demo`.

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

The result is the folder `dist/MSparser` (about 550 MB). The converters are
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
- Not generated: ions with a cation and an anion attached to the same
  molecule, a molecule combined with itself, ions of three or more molecules,
  neutral losses and radical ions.
- The ion list can be long, since every combination allowed by the rules is
  proposed; many of these ions will not exist in a given sample.
- Stability mode: with few traces, the limits derived from the data are not
  reliable; the program warns and the limits can be set by hand.
- Intensity mode: the program needs permission to write in the data folders;
  within a folder the same compound list applies to all files; files are
  processed one after another in the application.
- The stability mode still uses a separate project folder.
