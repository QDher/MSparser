"""Settings read from the data files, used when a setting is left on "auto".

Detected per acquisition at ingestion and stored in work/detected.json:
vendor, analyzer, polarity, scan range, scan rate, duration, resolving power.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

# used when the resolving power is neither in the file nor in the configuration;
# low values, so that close-lying targets are flagged rather than missed
DEFAULT_RESOLUTION = {
    "orbitrap": {"value": 60000.0, "at_mz": 200.0, "scaling": "orbitrap"},
    "tof": {"value": 20000.0, "at_mz": 200.0, "scaling": "tof"},
}
CANDIDATE_LENGTHS = [3, 5, 8, 10, 15, 20, 30, 45, 60, 90, 120]


def detect_vendor(scans: pd.DataFrame) -> str:
    sid = str(scans["native_id"].iloc[0]) if len(scans) else ""
    if "controllerType=" in sid:
        return "thermo"
    if "function=" in sid:
        return "waters"
    return "unknown"


def detect_file(scans: pd.DataFrame, cfg: dict) -> dict:
    """Summary of one acquisition. `scans` holds all scans with the 'kept' column."""
    vendor = cfg["vendor"] if cfg["vendor"] != "auto" else detect_vendor(scans)
    kept = scans[scans["kept"]].sort_values("time_s")
    filt = str(kept["filter_string"].dropna().iloc[0]) if kept["filter_string"].notna().any() else ""
    analyzer = "tof" if vendor == "waters" else "orbitrap"
    rp = kept["resolving_power"].dropna()
    t = kept["time_s"].to_numpy()
    return dict(
        vendor=vendor, analyzer=analyzer, filter_string=filt,
        polarity=int(kept["polarity"].iloc[0]) if len(kept) else 0,
        scan_range=[float(kept["win_lo"].iloc[0]), float(kept["win_hi"].iloc[0])] if len(kept) else [np.nan, np.nan],
        n_scans=int(len(kept)), n_dropped=int((~scans["kept"]).sum()),
        functions_kept=sorted(int(f) for f in kept["function"].unique()),
        duration_s=float(t[-1] - t[0]) if len(t) > 1 else 0.0,
        scan_interval_s=float(np.median(np.diff(t))) if len(t) > 1 else np.nan,
        resolving_power=float(rp.median()) if len(rp) else None,
    )


def save_detected(work: Path, detected: dict):
    (work / "detected.json").write_text(json.dumps(detected, indent=2))


def load_detected(work: Path) -> dict:
    p = work / "detected.json"
    return json.loads(p.read_text()) if p.exists() else {}


def resolve_resolution(cfg: dict, det: dict | None) -> dict:
    """Resolving power for one acquisition: configuration, else file, else default."""
    res = cfg["resolution"]
    det = det or {}
    base = dict(DEFAULT_RESOLUTION[det.get("analyzer", "orbitrap")])
    if isinstance(res, dict) and res.get("value") not in (None, "auto"):
        base.update({k: v for k, v in res.items() if v is not None}, source="configuration")
        if "scaling" not in res or res["scaling"] in (None, "auto"):
            base["scaling"] = DEFAULT_RESOLUTION[det.get("analyzer", "orbitrap")]["scaling"]
    elif det.get("resolving_power"):
        base.update(value=float(det["resolving_power"]), source="file")
    else:
        base["source"] = "default"
    return base


def resolve_windows(cfg: dict, detected: dict) -> dict:
    """Window lengths and step from scan rate and duration when set to "auto"."""
    w = dict(cfg["windows"])
    dets = [d for d in detected.values() if d.get("n_scans", 0) > 1]
    interval = float(np.nanmax([d["scan_interval_s"] for d in dets])) if dets else 1.0
    duration = float(np.nanmin([d["duration_s"] for d in dets])) if dets else 30.0
    source = {}
    if w["step_s"] in (None, "auto"):
        w["step_s"] = 1.0 if interval <= 1.0 else float(np.ceil(interval))
        source["step_s"] = "auto"
    if w["lengths_s"] in (None, "auto") or not w["lengths_s"]:
        shortest = w["min_scans"] * interval          # a window must hold min_scans scans
        longest = 0.6 * duration
        lengths = [float(x) for x in CANDIDATE_LENGTHS if shortest <= x <= longest]
        if len(lengths) < 2:  # very short or slow acquisitions: multiples of the minimum
            lengths = sorted({round(shortest * f, 1) for f in (1.0, 1.5, 2.0, 3.0) if shortest * f <= max(longest, shortest)})
        w["lengths_s"] = lengths
        source["lengths_s"] = "auto"
    w["source"] = source
    return w
