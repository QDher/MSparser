"""Configuration loading. One YAML file per vendor run."""
from __future__ import annotations

import copy
from pathlib import Path

import yaml

DEFAULTS = {
    "vendor": "thermo",  # thermo | waters
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
    "scan_filter": {
        "ms_level": 1,
        "functions": None,  # Waters: list of function numbers to keep, e.g. [1]
    },
    "extraction": {
        "ppm_windows": [10, 5],
        "min_window_mda": 0.5,          # absolute half-width floor at low m/z
        "min_consecutive_scans": 5,     # detection criterion
        "multi_centroid_fraction": 0.1,
        "isotope_ratio_tolerance": 0.5,  # relative deviation accepted
        "exclude_isotope_mismatch": True,
        "saturation_intensity": None,   # TOF: flag traces reaching this value
        "max_injection_time_ms": None,  # Orbitrap: configured maximum IT
    },
    "resolution": {"value": 60000, "at_mz": 200, "scaling": "orbitrap"},
    "windows": {
        "lengths_s": [3, 5, 8, 10, 15, 20],
        "step_s": 1.0,
        "min_scans": 8,
        "min_detected_fraction": 0.9,
    },
    "phases": {"plateau_fraction": 0.8, "onset_fraction": 0.1, "penalty_factor": 3.0},
    "thresholds": {
        "n_sigma": 3.0, "reference_length_s": None, "min_rsd_factor": 1.5,
        # optional manual values replacing the data-derived ones
        "override": {"rsd_factor": None, "slope_limit": None},
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
    if cfg["vendor"] not in ("thermo", "waters"):
        raise ValueError("vendor must be 'thermo' or 'waters'")
    cfg["extraction"]["ppm_windows"] = sorted(
        {float(p) for p in cfg["extraction"]["ppm_windows"]}, reverse=True
    )
    return cfg
