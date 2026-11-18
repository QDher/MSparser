"""Intensity-extraction mode.

For every target ion in every acquisition, report the signal intensity,
whatever its stability:

- the m/z window is adaptive: the configured widths are tried in increasing
  order and the first one containing a signal is used;
- inside that window the centroids are grouped into peaks; the peak closest to
  the theoretical m/z, the most intense one, and any peak of at least 10 % of
  the most intense are reported as candidates;
- the intensity of a peak is its mean over the most stable window of fixed
  length (lowest RSD among the sliding windows);
- low stability, variation between windows, etc. are flagged, never excluded.

Files are found recursively, so they can be organised in any folder tree.
Results are cached per file: a re-run only processes new or changed files.
"""
from __future__ import annotations

import hashlib
import json
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

from .describe import window_stats
from .detect import detect_file
from .extract import half_width
from .io import filter_scans, read_csv_spectra, read_mzml
from .targets import targets_for_file

KEY = ["file", "target_id"]
FLAGS = ["low_stability", "window_variation", "gaps", "multiple_peaks", "wide_mz_window",
         "isotope_mismatch", "saturated", "short_acquisition"]


# ------------------------------------------------------------------- files
def discover_files(data_dir: str | Path) -> pd.DataFrame:
    """mzML and CSV files below a folder, at any depth. `file` is the path
    relative to the folder without extension, so that names stay unique."""
    data_dir = Path(data_dir)
    rows = []
    if data_dir.is_dir():
        for p in sorted(data_dir.rglob("*")):
            if p.is_file() and p.suffix.lower() in (".mzml", ".csv"):
                rel = p.relative_to(data_dir).with_suffix("")
                rows.append(dict(file=rel.as_posix(), folder=rel.parent.as_posix() if str(rel.parent) != "." else "",
                                 name=rel.name, path=str(p)))
    return pd.DataFrame(rows, columns=["file", "folder", "name", "path"])


# ------------------------------------------------------------------- peaks
def group_peaks(mz: np.ndarray, center: float, width_ppm: float) -> np.ndarray:
    """Group centroids collected over many scans into peaks. Returns a label per
    centroid. Peaks are the modes of the m/z density (Gaussian kernel of
    `width_ppm`); each centroid is assigned to the nearest mode."""
    if len(mz) == 0:
        return np.array([], dtype=int)
    ppm = (mz - center) / center * 1e6
    lo, hi = ppm.min() - 3 * width_ppm, ppm.max() + 3 * width_ppm
    grid = np.arange(lo, hi + 0.25 * width_ppm, 0.25 * width_ppm)
    dens = np.exp(-0.5 * ((grid[:, None] - ppm[None, :]) / width_ppm) ** 2).sum(axis=1)
    is_max = np.r_[True, dens[1:] > dens[:-1]] & np.r_[dens[:-1] >= dens[1:], True]
    modes = grid[is_max & (dens >= 0.02 * dens.max())]
    if len(modes) == 0:
        modes = np.array([np.median(ppm)])
    return np.argmin(np.abs(ppm[:, None] - modes[None, :]), axis=1)


def _longest_run(mask: np.ndarray) -> int:
    best = run = 0
    for m in mask:
        run = run + 1 if m else 0
        best = max(best, run)
    return best


def _candidates(t, slices, p_mz, p_int, target_mz, tol, ic):
    """Peaks in the m/z window: list of (trace intensity, trace m/z) per peak."""
    idx_scan, idx_peak = [], []
    for i, (a, b) in enumerate(slices):
        lo = a + np.searchsorted(p_mz[a:b], target_mz - tol, "left")
        hi = a + np.searchsorted(p_mz[a:b], target_mz + tol, "right")
        idx_scan.extend([i] * (hi - lo))
        idx_peak.extend(range(lo, hi))
    if not idx_peak:
        return []
    idx_scan, idx_peak = np.array(idx_scan), np.array(idx_peak)
    mz, inten = p_mz[idx_peak], p_int[idx_peak]
    labels = group_peaks(mz, target_mz, ic["peak_width_ppm"])
    out = []
    for lab in np.unique(labels):
        m = labels == lab
        center = np.median(mz[m])
        y, ymz = np.full(len(t), np.nan), np.full(len(t), np.nan)
        # one value per scan: if a scan has two centroids in the peak, keep the closest to its centre
        order = np.argsort(np.abs(mz[m] - center))[::-1]
        y[idx_scan[m][order]] = inten[m][order]
        ymz[idx_scan[m][order]] = mz[m][order]
        if _longest_run(np.isfinite(y)) >= ic["min_scans"]:
            out.append((y, ymz))
    return out


def _measure(t, y, ic) -> tuple[dict, pd.DataFrame]:
    """Intensity of one peak: mean over the most stable window of fixed length.
    Returns (values, table of all windows)."""
    length, step = float(ic["window_length_s"]), float(ic["window_step_s"])
    w = {"min_scans": max(3, int(ic["min_scans"])), "min_detected_fraction": 0.0}
    rows = []
    k = 0
    while t[0] + k * step + length <= t[-1] + 1e-9:
        s = t[0] + k * step
        r = window_stats(t, y, s, s + length, w)
        if r["n"] >= w["min_scans"] and np.isfinite(r["rsd"]):
            rows.append(dict(start=round(s, 3), end=round(s + length, 3), **r))
        k += 1
    wins = pd.DataFrame(rows)
    short = wins.empty
    if short:  # acquisition shorter than the window, or too few scans: use everything
        r = window_stats(t, y, t[0], t[-1] + 1e-6, w)
        wins = pd.DataFrame([dict(start=round(float(t[0]), 3), end=round(float(t[-1]), 3), **r)])
    best = wins.loc[wins["rsd"].idxmin()] if wins["rsd"].notna().any() else wins.iloc[0]
    means = wins["mean"].dropna()
    variation = float((means.max() - means.min()) / best["mean"] * 100) if len(means) > 1 and best["mean"] > 0 else 0.0
    det = np.isfinite(y)
    vals = dict(
        intensity=float(best["mean"]), rsd_pct=float(best["rsd"] * 100), slope_pct_per_s=float(best["rel_slope"]),
        window_start_s=float(best["start"]), window_end_s=float(best["end"]),
        scans_in_window=int(best["n"]), detected_fraction_in_window=float(best["n"] / max(best["n_scans"], 1)),
        window_variation_pct=variation, window_min_intensity=float(means.min()) if len(means) else np.nan,
        window_max_intensity=float(means.max()) if len(means) else np.nan, n_windows=int(len(means)),
        scans_detected=int(det.sum()), mean_all_scans=float(np.nanmean(y)), short_acquisition=bool(short),
    )
    wins = wins.assign(selected=wins.index == best.name)[["start", "end", "n", "mean", "rsd", "rel_slope", "selected"]]
    return vals, wins


# ------------------------------------------------------------ one acquisition
def process_file(path: str, file: str, targets: pd.DataFrame, cfg: dict) -> dict:
    """All candidate peaks of all target ions in one acquisition."""
    ic, ex = cfg["intensity"], cfg["extraction"]
    p = Path(path)
    scans, peaks = read_csv_spectra(p, cfg) if p.suffix.lower() == ".csv" else read_mzml(p)
    scans = filter_scans(scans, cfg)
    info = detect_file(scans, cfg)
    kept = scans[scans["kept"]].sort_values("time_s").reset_index(drop=True)
    empty = dict(info=info, candidates=pd.DataFrame(), windows=pd.DataFrame(), traces=pd.DataFrame())
    if len(kept) < 2:
        return empty
    peaks = peaks[peaks["scan_idx"].isin(kept["scan_idx"])].sort_values(["scan_idx", "mz"])
    p_idx, p_mz = peaks["scan_idx"].to_numpy(), peaks["mz"].to_numpy()
    p_int = peaks["intensity"].to_numpy(dtype=float)
    slices = [(np.searchsorted(p_idx, s, "left"), np.searchsorted(p_idx, s, "right")) for s in kept["scan_idx"]]
    t = kept["time_s"].to_numpy()
    win_lo, win_hi = kept["win_lo"].iloc[0], kept["win_hi"].iloc[0]
    steps = sorted(float(x) for x in ic["ppm_steps"])
    cands, windows, traces = [], [], []
    for tg in targets_for_file(targets, None, info["polarity"]).itertuples(index=False):
        base = dict(file=file, target_id=tg.target_id, component=tg.component, ion=tg.ion,
                    type=getattr(tg, "type", ""), charge=tg.charge, mz_theoretical=tg.mz)
        if np.isfinite(win_lo) and np.isfinite(win_hi) and not win_lo <= tg.mz <= win_hi:
            cands.append({**base, "status": "out of scan range"})
            continue
        found, ppm_used = [], np.nan
        for ppm in steps:  # adaptive m/z window
            tol = half_width(tg.mz, ppm, ex["min_window_mda"])
            found = _candidates(t, slices, p_mz, p_int, tg.mz, tol, ic)
            if found:
                ppm_used = ppm
                break
        if not found:
            cands.append({**base, "status": "not detected", "mz_window_ppm": steps[-1]})
            continue
        measured = []
        for y, ymz in found:
            vals, wins = _measure(t, y, ic)
            mz_obs = float(np.nanmedian(ymz))
            iso_ratio = np.nan
            if np.isfinite(tg.iso_ratio):  # isotope peak, at the same offset from the observed m/z
                shift = tg.iso_mz - tg.mz
                itol = half_width(tg.mz, 3 * ic["peak_width_ppm"], ex["min_window_mda"])
                iso = np.full(len(t), np.nan)
                for i, (a, b) in enumerate(slices):
                    lo = a + np.searchsorted(p_mz[a:b], mz_obs + shift - itol, "left")
                    hi = a + np.searchsorted(p_mz[a:b], mz_obs + shift + itol, "right")
                    if hi > lo:
                        iso[i] = p_int[lo:hi].max()
                both = np.isfinite(iso) & np.isfinite(y)
                if both.sum() >= 3:
                    iso_ratio = float(np.median(iso[both] / y[both]))
            measured.append((vals, wins, y, ymz, mz_obs, iso_ratio))
        top = max(m[0]["intensity"] for m in measured)
        errs = [abs(m[4] - tg.mz) for m in measured]
        i_close, i_top = int(np.argmin(errs)), int(np.argmax([m[0]["intensity"] for m in measured]))
        keep = [i for i, m in enumerate(measured)
                if i in (i_close, i_top) or m[0]["intensity"] >= ic["min_relative_intensity"] * top]
        for i in keep:
            vals, wins, y, ymz, mz_obs, iso_ratio = measured[i]
            sat = ex["saturation_intensity"]
            iso_dev = abs(iso_ratio / tg.iso_ratio - 1) if np.isfinite(iso_ratio) and tg.iso_ratio > 0 else np.nan
            row = {
                **base, "status": "detected", "mz_window_ppm": ppm_used, "mz_observed": mz_obs,
                "ppm_error": (mz_obs - tg.mz) / tg.mz * 1e6, **vals,
                "relative_to_most_intense_pct": vals["intensity"] / top * 100,
                "is_closest": i == i_close, "is_most_intense": i == i_top, "n_peaks_in_window": len(keep),
                "isotope_ratio_theory": tg.iso_ratio, "isotope_ratio_observed": iso_ratio,
                "low_stability": vals["rsd_pct"] > ic["rsd_flag_pct"],
                "window_variation": vals["window_variation_pct"] >= ic["variation_flag_pct"],
                "gaps": vals["detected_fraction_in_window"] < 0.9,
                "multiple_peaks": len(keep) > 1,
                "wide_mz_window": ppm_used > steps[0],
                "isotope_mismatch": bool(np.isfinite(iso_dev) and iso_dev > ex["isotope_ratio_tolerance"]),
                "saturated": bool(sat is not None and np.nanmax(y) >= sat),
            }
            cands.append(row)
            ident = dict(file=file, target_id=tg.target_id, mz_observed=mz_obs)
            windows.append(wins.assign(**ident))
            d = np.isfinite(y)
            traces.append(pd.DataFrame(dict(time_s=t[d], intensity=y[d], mz=ymz[d])).assign(**ident))
    cat = lambda parts: pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()
    return dict(info=info, candidates=pd.DataFrame(cands), windows=cat(windows), traces=cat(traces))


def _cache_key(path: Path, targets: pd.DataFrame, cfg: dict) -> str:
    st = path.stat()
    t = targets[targets["use"] == 1][["target_id", "mz", "iso_mz", "iso_ratio"]].round(6).to_json()
    # settings that do not change what is measured in a file are left out of the key
    ic = {k: v for k, v in cfg["intensity"].items() if k not in ("default_peak", "jobs")}
    params = json.dumps([ic, cfg["extraction"], cfg["scan_filter"], cfg["csv"], cfg["vendor"]],
                        sort_keys=True, default=str)
    return hashlib.sha1(f"{path}|{st.st_size}|{st.st_mtime_ns}|{t}|{params}".encode()).hexdigest()[:20]


def _one(args):
    path, file, targets, cfg, cache = args
    key = _cache_key(Path(path), targets, cfg)
    base = Path(cache) / key
    if (base.with_suffix(".json")).exists():
        return file, "cached", None, key
    try:
        r = process_file(path, file, targets, cfg)
        for name in ("candidates", "windows", "traces"):
            r[name].to_parquet(f"{base}.{name}.parquet", index=False)
        base.with_suffix(".json").write_text(json.dumps(r["info"]))
        return file, "processed", None, key
    except Exception as e:  # one unreadable file must not stop the batch
        return file, "failed", f"{type(e).__name__}: {e}", key


# ---------------------------------------------------------------- the batch
def run_batch(cfg: dict, targets: pd.DataFrame, progress=None) -> dict:
    """Process every file below the data folder (cached per file) and gather
    the results. `progress(done, total, file)` is called after each file."""
    work = Path(cfg["paths"]["work_dir"])
    cache = work / "intensity_cache"
    cache.mkdir(parents=True, exist_ok=True)
    files = discover_files(cfg["paths"]["mzml_dir"])
    if files.empty:
        raise FileNotFoundError(f"no mzML or CSV file found in {cfg['paths']['mzml_dir']} or its subfolders")
    jobs = [(r.path, r.file, targets, cfg, str(cache)) for r in files.itertuples(index=False)]
    n_jobs = max(1, int(cfg["intensity"].get("jobs") or 1))
    results = []
    if n_jobs > 1:
        with ProcessPoolExecutor(max_workers=n_jobs) as pool:
            for i, res in enumerate(pool.map(_one, jobs), 1):
                results.append(res)
                if progress:
                    progress(i, len(jobs), res[0])
    else:
        for i, job in enumerate(jobs, 1):
            res = _one(job)
            results.append(res)
            if progress:
                progress(i, len(jobs), res[0])
    status = pd.DataFrame(results, columns=["file", "processing", "error", "key"]).merge(files, on="file")
    parts = {"candidates": [], "windows": [], "traces": []}
    infos = []
    for r in status.itertuples(index=False):
        if r.processing == "failed":
            infos.append(dict(file=r.file))
            continue
        base = cache / r.key
        infos.append(dict(file=r.file, **json.loads(base.with_suffix(".json").read_text())))
        for name in parts:
            df = pd.read_parquet(f"{base}.{name}.parquet")
            if len(df):
                parts[name].append(df)
    cat = lambda p: pd.concat(p, ignore_index=True) if p else pd.DataFrame()
    info = pd.DataFrame(infos)
    if "scan_range" in info:
        info["scan_range"] = info["scan_range"].map(lambda v: "-".join(f"{x:g}" for x in v) if isinstance(v, list) else "")
    if "functions_kept" in info:
        info["functions_kept"] = info["functions_kept"].map(lambda v: ",".join(map(str, v)) if isinstance(v, list) else "")
    files_out = status.drop(columns=["key", "path"]).merge(info, on="file", how="left")
    return dict(files=files_out, candidates=cat(parts["candidates"]), windows=cat(parts["windows"]),
                traces=cat(parts["traces"]))


# ------------------------------------------------------- choice of the peak
def apply_selection(cands: pd.DataFrame, selection: pd.DataFrame | None, cfg: dict) -> pd.DataFrame:
    """Mark the peak used for each ion and file (`selected`, `selected_by`).
    `selection` holds the user's choices (file, target_id, mz_observed); the
    other ions follow the default rule."""
    c = cands.copy()
    c["selected"], c["selected_by"] = False, ""
    det = c["status"] == "detected"
    rule = "is_most_intense" if cfg["intensity"]["default_peak"] == "most_intense" else "is_closest"
    label = "rule: most intense" if rule == "is_most_intense" else "rule: closest m/z"
    chosen = {}
    if selection is not None and len(selection):
        for r in selection.itertuples(index=False):
            chosen[(r.file, r.target_id)] = float(r.mz_observed)
    for key, g in c[det].groupby(KEY, sort=False):
        idx = None
        if key in chosen:
            tol = 2 * cfg["intensity"]["peak_width_ppm"] * 1e-6 * g["mz_theoretical"].iloc[0]
            d = (g["mz_observed"] - chosen[key]).abs()
            if d.min() <= tol:
                idx, by = d.idxmin(), "user"
        if idx is None:
            idx, by = g.index[g[rule].astype(bool)][0], label
        c.loc[idx, ["selected", "selected_by"]] = [True, by]
    return c


def build_outputs(cands: pd.DataFrame, files: pd.DataFrame) -> dict:
    """Long table (one row per file and ion) and wide table (files x ions)."""
    if cands.empty:
        return dict(long=pd.DataFrame(), matrix=pd.DataFrame())
    sel = cands[cands["selected"] | (cands["status"] != "detected")].copy()
    flag_cols = [f for f in FLAGS if f in sel]
    sel["flags"] = sel[flag_cols].fillna(False).astype(bool).apply(
        lambda r: "; ".join(f for f in flag_cols if r[f]), axis=1)
    sel = sel.merge(files[["file", "folder", "name"]], on="file", how="left")
    first = ["folder", "name", "file", "target_id", "status", "intensity", "flags", "rsd_pct", "window_variation_pct",
             "mz_window_ppm", "ppm_error", "mz_theoretical", "mz_observed", "window_start_s", "window_end_s"]
    long = sel[[c for c in first if c in sel] + [c for c in sel.columns if c not in first]]
    long = long.sort_values(["file", "target_id"]).reset_index(drop=True)
    det = long[long["status"] == "detected"]
    matrix = det.pivot(index="file", columns="target_id", values="intensity")
    matrix = matrix.reindex(sorted(files["file"])).reset_index()
    matrix = files[["file", "folder", "name"]].merge(matrix, on="file", how="right")
    return dict(long=long, matrix=matrix)
