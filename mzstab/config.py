"""Configuration loading. One YAML file per vendor run."""
from __future__ import annotations

import copy
from pathlib import Path

import yaml

DEFAULTS = {
    # Settings on "auto" are read from the data files (see detect.py).
    "vendor": "auto",  # auto | thermo | waters
    "paths": {
        "samples": "samples.csv",
        "components": "components.yaml",
        "raw_dir": "raw",
        "mzml_dir": "mzml",
        "work_dir": "work",
        "report_dir": "report",
    },
    "convert": {
        "thermo_parser": "ThermoRawFileParser",
        "msconvert": "msconvert",
        "lockmass_mz": None,       # Waters: lock-mass m/z (positive mode)
        "lockmass_mz_neg": None,   # Waters: lock-mass m/z (negative mode)
        "lockmass_tol_da": 0.5,
    },
    # CSV input: one row per centroid peak, columns time, m/z, intensity and
    # optionally scan. Column names are detected; set them here if detection fails.
    "csv": {
        "time_unit": "s",       # s | min
        "polarity": 1,          # +1 or -1 (not stored in CSV)
        "time_col": None, "mz_col": None, "intensity_col": None, "scan_col": None,
    },
    "targets": {"pairs": True},  # also propose ions made of two different components
    "scan_filter": {
        "ms_level": 1,
        "functions": None,  # Waters functions to keep; None = function 1 for Waters files
    },
    "extraction": {
        "ppm_windows": [15],            # several values are all run and compared
        "min_window_mda": 0.5,          # absolute half-width floor at low m/z
        "min_consecutive_scans": 5,     # detection criterion
        "multi_centroid_fraction": 0.1,
        "isotope_ratio_tolerance": 0.5,  # relative deviation accepted
        "exclude_isotope_mismatch": True,
        "saturation_intensity": None,   # TOF: flag traces reaching this value
        "max_injection_time_ms": None,  # Orbitrap: configured maximum IT
    },
    "resolution": {"value": "auto", "at_mz": 200},  # auto: from the file, else a default
    "windows": {
        "lengths_s": "auto",  # auto: from scan rate and acquisition duration
        "step_s": "auto",
        "min_scans": 8,
        "min_detected_fraction": 0.9,
    },
    "phases": {"plateau_fraction": 0.8, "onset_fraction": 0.1, "penalty_factor": 3.0},
    "thresholds": {
        "n_sigma": 3.0, "reference_length_s": None, "min_rsd_factor": 1.5,
        # Stability limits. "auto" derives the value from the data (spread of the
        # most stable window of each trace); a number is used as given.
        "rsd_factor": "auto",   # RSD limit = this factor x the noise model
        "slope_limit": "auto",  # maximum |slope|, in % of the mean intensity per second
    },
    # Windows set by the user, on acquisition time, as [start_s, end_s].
    # "targets" maps a target_id to its own window and takes precedence.
    "manual_windows": {"global": None, "targets": {}},
    # Intensity-extraction mode (see intensity.py)
    "mode": "stability",  # stability | intensity
    "intensity": {
        "ppm_steps": [5, 10, 15, 20],   # tried in this order; the first with a signal is used
        "min_scans": 5,                 # consecutive scans required to report a peak
        "window_length_s": 10.0,        # fixed window; the most stable one is used
        "window_step_s": 1.0,
        "variation_flag_pct": 25.0,     # flag when 10 s windows differ by this much
        "rsd_flag_pct": 30.0,           # flag "low stability" above this RSD
        "peak_width_ppm": 1.0,          # m/z scatter of one peak between scans
        "min_relative_intensity": 0.10, # other peaks reported if >= this x the most intense
        "default_peak": "closest",      # closest | most_intense
        "jobs": 1,                      # parallel processes (command line only)
    },
    "validation": {"test_fraction": 0.5, "seed": 0, "min_files": 4},
}


def _merge(base: dict, over: dict) -> dict:
    for k, v in over.items():
        if isinstance(v, dict) and isinstance(base.get(k), dict):
            _merge(base[k], v)
        else:
            base[k] = v
    return base


def load_config(path: str | Path) -> dict:
    path = Path(path).resolve()
    with open(path) as fh:
        user = yaml.safe_load(fh) or {}
    cfg = _merge(copy.deepcopy(DEFAULTS), user)
    root = path.parent
    cfg["paths"] = {k: str((root / v).resolve()) for k, v in cfg["paths"].items()}
    if cfg["vendor"] not in ("auto", "thermo", "waters"):
        raise ValueError("vendor must be 'auto', 'thermo' or 'waters'")
    cfg["extraction"]["ppm_windows"] = sorted(
        {float(p) for p in cfg["extraction"]["ppm_windows"]}, reverse=True
    )
    return cfg
