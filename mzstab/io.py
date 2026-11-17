"""mzML reading, scan filtering, Parquet storage."""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd
from pyteomics import mzml

from .detect import detect_file, detect_vendor


def _seconds(value) -> float:
    unit = getattr(value, "unit_info", None)
    return float(value) * 60.0 if unit == "minute" else float(value)


def read_mzml(path: str | Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return (scans, peaks). peaks: scan_idx, mz, intensity (zeros removed)."""
    scans, mz_parts, int_parts, idx_parts = [], [], [], []
    with mzml.read(str(path)) as reader:
        for i, s in enumerate(reader):
            scan = (s.get("scanList", {}).get("scan") or [{}])[0]
            win = (scan.get("scanWindowList", {}).get("scanWindow") or [{}])[0]
            sid = s.get("id", "")
            fn = re.search(r"function=(\d+)", sid)
            mz = np.asarray(s.get("m/z array", []), dtype=np.float64)
            inten = np.asarray(s.get("intensity array", []), dtype=np.float64)
            keep = inten > 0
            mz, inten = mz[keep], inten[keep]
            order = np.argsort(mz, kind="stable")
            mz, inten = mz[order], inten[order]
            it = scan.get("ion injection time")
            scans.append(
                dict(
                    scan_idx=i, native_id=sid,
                    function=int(fn.group(1)) if fn else 1,
                    ms_level=int(s.get("ms level", 0)),
                    polarity=1 if "positive scan" in s else (-1 if "negative scan" in s else 0),
                    time_s=_seconds(scan["scan start time"]) if "scan start time" in scan else np.nan,
                    tic=float(s.get("total ion current", inten.sum())),
                    injection_time_ms=float(it) if it is not None else np.nan,
                    win_lo=float(win.get("scan window lower limit", np.nan)),
                    win_hi=float(win.get("scan window upper limit", np.nan)),
                    centroided="profile spectrum" not in s,
                    filter_string=scan.get("filter string"),
                    resolving_power=float(scan["mass resolving power"]) if "mass resolving power" in scan else np.nan,
                    n_peaks=len(mz),
                )
            )
            mz_parts.append(mz)
            int_parts.append(inten.astype(np.float32))
            idx_parts.append(np.full(len(mz), i, dtype=np.int32))
    scans = pd.DataFrame(scans)
    peaks = pd.DataFrame(
        dict(
            scan_idx=np.concatenate(idx_parts) if idx_parts else np.array([], np.int32),
            mz=np.concatenate(mz_parts) if mz_parts else np.array([]),
            intensity=np.concatenate(int_parts) if int_parts else np.array([], np.float32),
        )
    )
    return scans, peaks


CSV_ALIASES = {
    "time_col": ("time", "time_s", "rt", "retention time", "retention_time", "scan time",
                 "scan_time", "time (s)", "time (min)", "rt (min)", "rt (s)"),
    "mz_col": ("mz", "m/z", "mass", "mass/charge"),
    "intensity_col": ("intensity", "int", "abundance", "counts", "signal"),
    "scan_col": ("scan", "scan_idx", "scan number", "scan_number", "spectrum"),
}


def read_csv_spectra(path: str | Path, cfg: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Centroided spectra from a long-format CSV: one row per peak with time,
    m/z and intensity (and optionally a scan number). Rows are grouped into
    scans by the scan column if present, otherwise by identical time values."""
    c = cfg["csv"]
    df = pd.read_csv(path, sep=None, engine="python")
    lower = {str(col).strip().lower(): col for col in df.columns}
    cols = {}
    for key, aliases in CSV_ALIASES.items():
        if c.get(key):
            if c[key] not in df.columns:
                raise ValueError(f"{Path(path).name}: column '{c[key]}' not found")
            cols[key] = c[key]
        else:
            cols[key] = next((lower[a] for a in aliases if a in lower), None)
    missing = [k for k in ("time_col", "mz_col", "intensity_col") if cols[k] is None]
    if missing:
        raise ValueError(
            f"{Path(path).name}: cannot identify {', '.join(missing)} among columns "
            f"{list(df.columns)}; set csv.{missing[0]} in the configuration"
        )
    d = pd.DataFrame({
        "time_s": pd.to_numeric(df[cols["time_col"]], errors="coerce") * (60.0 if c["time_unit"] == "min" else 1.0),
        "mz": pd.to_numeric(df[cols["mz_col"]], errors="coerce"),
        "intensity": pd.to_numeric(df[cols["intensity_col"]], errors="coerce"),
    })
    key = df[cols["scan_col"]] if cols["scan_col"] else d["time_s"]
    d["key"] = key.to_numpy()
    d = d.dropna(subset=["time_s", "mz", "intensity", "key"])
    first_time = d.groupby("key")["time_s"].min().sort_values(kind="stable")
    d["scan_idx"] = d["key"].map({k: i for i, k in enumerate(first_time.index)}).astype(np.int32)
    pk = d[d["intensity"] > 0].sort_values(["scan_idx", "mz"])
    tic = d.groupby("scan_idx")["intensity"].sum()
    n_pk = pk.groupby("scan_idx").size()
    idx = np.arange(len(first_time))
    scans = pd.DataFrame(dict(
        scan_idx=idx, native_id=[f"scan={i + 1}" for i in idx], function=1, ms_level=1,
        polarity=int(c["polarity"]), time_s=first_time.to_numpy(),
        tic=tic.reindex(idx).fillna(0.0).to_numpy(), injection_time_ms=np.nan,
        win_lo=np.nan, win_hi=np.nan, centroided=True, filter_string=None, resolving_power=np.nan,
        n_peaks=n_pk.reindex(idx).fillna(0).astype(int).to_numpy(),
    ))
    peaks = pd.DataFrame(dict(scan_idx=pk["scan_idx"].to_numpy(), mz=pk["mz"].to_numpy(),
                              intensity=pk["intensity"].to_numpy(dtype=np.float32)))
    return scans, peaks


def find_input(folder: str | Path, stem: str) -> Path:
    """The spectra file of an acquisition: <stem>.mzML, or <stem>.csv."""
    for ext in (".mzML", ".mzml", ".csv", ".CSV"):
        p = Path(folder) / f"{stem}{ext}"
        if p.exists():
            return p
    raise FileNotFoundError(f"no {stem}.mzML or {stem}.csv in {folder}")


def filter_scans(scans: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    """Add 'kept' and 'drop_reason'. Keeps the configured MS level (and Waters
    functions), then the dominant polarity and the dominant scan range."""
    sf = cfg["scan_filter"]
    scans = scans.copy()
    reason = pd.Series("", index=scans.index)
    reason[scans["ms_level"] != sf["ms_level"]] = "ms_level"
    functions = sf.get("functions")
    if not functions and detect_vendor(scans) == "waters" and cfg["vendor"] != "thermo":
        functions = [1]  # Waters: function 1 is the MS function; others are lock mass or high energy
    if functions:
        reason[(reason == "") & ~scans["function"].isin(functions)] = "function"
    ok = reason == ""
    if ok.any():
        pol = scans.loc[ok, "polarity"].mode().iloc[0]
        reason[ok & (scans["polarity"] != pol)] = "polarity"
        ok = reason == ""
        w = scans.loc[ok, ["win_lo", "win_hi"]].round(1).fillna(-1.0)
        rng = w["win_lo"].map(str) + "-" + w["win_hi"].map(str)
        dom = rng.mode().iloc[0]
        reason[rng.index[rng != dom]] = "scan_range"
    scans["drop_reason"] = reason
    scans["kept"] = reason == ""
    return scans


def ingest_file(path: Path, cfg: dict, work: Path) -> tuple[pd.DataFrame, dict]:
    is_csv = path.suffix.lower() == ".csv"
    scans, peaks = read_csv_spectra(path, cfg) if is_csv else read_mzml(path)
    scans = filter_scans(scans, cfg)
    if not scans["centroided"].all():
        print(f"WARNING {path.name}: profile spectra found; convert with peak picking")
    peaks = peaks[peaks["scan_idx"].isin(scans.loc[scans["kept"], "scan_idx"])]
    scans.insert(0, "file", path.stem)
    (work / "scans").mkdir(parents=True, exist_ok=True)
    (work / "peaks").mkdir(parents=True, exist_ok=True)
    scans.to_parquet(work / "scans" / f"{path.stem}.parquet", index=False)
    peaks.to_parquet(work / "peaks" / f"{path.stem}.parquet", index=False)
    return scans, detect_file(scans, cfg)


def load_samples(cfg: dict) -> pd.DataFrame:
    """samples.csv: file, reaction, components (';'-separated, empty = all), run_order."""
    df = pd.read_csv(cfg["paths"]["samples"], dtype=str, keep_default_na=False)
    if "file" not in df:
        raise ValueError("samples.csv needs a 'file' column")
    df["file"] = df["file"].map(lambda f: Path(f).stem)
    if df["file"].duplicated().any():
        raise ValueError("duplicate file names in samples.csv")
    for col in ("reaction", "components"):
        if col not in df:
            df[col] = ""
    order = df["run_order"] if "run_order" in df else pd.Series("", index=df.index)
    df["run_order"] = pd.to_numeric(order, errors="coerce")
    df["run_order"] = df["run_order"].fillna(pd.Series(np.arange(1, len(df) + 1), index=df.index))
    return df
