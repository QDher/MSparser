"""Targeted extraction of ion traces, QC flags, m/z window comparison."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from .targets import targets_for_file

MAD = 1.4826


def half_width(mz: float, ppm: float, min_mda: float) -> float:
    return max(mz * ppm * 1e-6, min_mda * 1e-3)


def resolving_power(mz: float, res: dict) -> float:
    if res["scaling"] == "orbitrap":
        return res["value"] * np.sqrt(res["at_mz"] / mz)
    return res["value"]


def _closest(mz_arr, int_arr, target, tol):
    lo = np.searchsorted(mz_arr, target - tol, "left")
    hi = np.searchsorted(mz_arr, target + tol, "right")
    if hi <= lo:
        return np.nan, np.nan, 0
    j = lo + int(np.argmin(np.abs(mz_arr[lo:hi] - target)))
    return float(int_arr[j]), float(mz_arr[j]), hi - lo


def _longest_run(mask: np.ndarray) -> int:
    best = run = 0
    for m in mask:
        run = run + 1 if m else 0
        best = max(best, run)
    return best


def extract_file(file: str, cfg: dict, targets: pd.DataFrame, injected, work: Path):
    """Return (traces, summary) for one acquisition, all ppm settings."""
    ex = cfg["extraction"]
    scans = pd.read_parquet(work / "scans" / f"{file}.parquet")
    scans = scans[scans["kept"]].sort_values("time_s").reset_index(drop=True)
    peaks = pd.read_parquet(work / "peaks" / f"{file}.parquet").sort_values(["scan_idx", "mz"])
    if scans.empty:
        return pd.DataFrame(), pd.DataFrame()
    polarity = int(scans["polarity"].iloc[0])
    tg = targets_for_file(targets, injected, polarity)

    p_idx = peaks["scan_idx"].to_numpy()
    p_mz = peaks["mz"].to_numpy()
    p_int = peaks["intensity"].to_numpy(dtype=float)
    bounds = [
        (np.searchsorted(p_idx, s, "left"), np.searchsorted(p_idx, s, "right"))
        for s in scans["scan_idx"]
    ]
    noise = float(np.percentile(p_int, 5)) if len(p_int) else np.nan
    win_lo, win_hi = scans["win_lo"].iloc[0], scans["win_hi"].iloc[0]
    it = scans["injection_time_ms"].to_numpy()
    it_at_max = False
    if np.isfinite(it).any():
        # AGC target not reached: injection time pinned at its maximum. Without
        # a configured maximum, the largest observed value is used as a proxy.
        max_it = ex["max_injection_time_ms"] or np.nanmax(it)
        it_at_max = bool(np.mean(it >= 0.99 * max_it) > 0.5)

    in_range = tg[
        ~((tg["mz"] < win_lo) | (tg["mz"] > win_hi))
    ] if np.isfinite(win_lo) and np.isfinite(win_hi) else tg
    unresolved = set()
    mzs = in_range["mz"].to_numpy()
    ids = in_range["target_id"].to_numpy()
    for a in range(len(mzs)):
        for b in range(a + 1, len(mzs)):
            if abs(mzs[a] - mzs[b]) < mzs[a] / resolving_power(mzs[a], cfg["resolution"]):
                unresolved.update((ids[a], ids[b]))

    traces, summary = [], []
    for t in tg.itertuples(index=False):
        out_of_range = t.target_id not in set(ids)
        for ppm in ex["ppm_windows"]:
            base = dict(file=file, target_id=t.target_id, component=t.component, ion=t.ion,
                        charge=t.charge, mz=t.mz, ppm=ppm, polarity=polarity,
                        n_scans=len(scans), out_of_range=out_of_range)
            if out_of_range:
                summary.append({**base, "not_detected": True, "analyzed": False})
                continue
            tol = half_width(t.mz, ppm, ex["min_window_mda"])
            rec = np.array([_closest(p_mz[a:b], p_int[a:b], t.mz, tol) for a, b in bounds])
            inten, mz_obs, n_cent = rec[:, 0], rec[:, 1], rec[:, 2].astype(int)
            iso = np.full(len(scans), np.nan)
            if np.isfinite(t.iso_ratio):
                iso = np.array([_closest(p_mz[a:b], p_int[a:b], t.iso_mz, tol)[0] for a, b in bounds])
            err = (mz_obs - t.mz) / t.mz * 1e6
            det = np.isfinite(inten)
            traces.append(pd.DataFrame(dict(
                file=file, target_id=t.target_id, ppm=ppm, scan_idx=scans["scan_idx"].to_numpy(),
                time_s=scans["time_s"].to_numpy(), intensity=inten, mz_obs=mz_obs,
                ppm_error=err, n_centroids=n_cent, iso_intensity=iso,
            )))
            detected = _longest_run(det) >= ex["min_consecutive_scans"]
            med = float(np.nanmedian(inten)) if det.any() else np.nan
            # isotope pattern check
            both = det & np.isfinite(iso)
            ratio_obs, iso_status = np.nan, "not_applicable"
            if np.isfinite(t.iso_ratio) and detected:
                if both.sum() >= 3:
                    ratio_obs = float(np.median(iso[both] / inten[both]))
                    dev = abs(ratio_obs / t.iso_ratio - 1)
                    iso_status = "ok" if dev <= ex["isotope_ratio_tolerance"] else "mismatch"
                elif med * t.iso_ratio < 3 * noise:
                    iso_status = "not_checkable"
                else:
                    iso_status = "mismatch"
            sat = ex["saturation_intensity"]
            mismatch = iso_status == "mismatch"
            summary.append({
                **base,
                "n_detected": int(det.sum()),
                "detected_fraction": float(det.mean()),
                "longest_run": _longest_run(det),
                "median_intensity": med,
                "max_intensity": float(np.nanmax(inten)) if det.any() else np.nan,
                "ppm_error_median": float(np.nanmedian(err)) if det.any() else np.nan,
                "ppm_error_sd": float(MAD * np.nanmedian(np.abs(err - np.nanmedian(err)))) if det.any() else np.nan,
                "multi_fraction": float(np.mean(n_cent > 1)),
                "iso_ratio_theory": t.iso_ratio,
                "iso_ratio_observed": ratio_obs,
                "iso_status": iso_status,
                "not_detected": not detected,
                "saturated": bool(sat is not None and det.any() and np.nanmax(inten) >= sat),
                "it_at_max": it_at_max,
                "unresolved": t.target_id in unresolved,
                "multiple_centroids": bool(np.mean(n_cent > 1) > ex["multi_centroid_fraction"]),
                "isotope_mismatch": mismatch,
                "analyzed": detected and not (mismatch and ex["exclude_isotope_mismatch"]),
            })
    traces = pd.concat(traces, ignore_index=True) if traces else pd.DataFrame()
    return traces, pd.DataFrame(summary)


def compare_windows(traces: pd.DataFrame, summary: pd.DataFrame) -> pd.DataFrame:
    """Per trace: effect of the narrow m/z window relative to the wide one."""
    ppms = sorted(traces["ppm"].unique(), reverse=True)
    if len(ppms) < 2:
        return pd.DataFrame()
    wide, narrow = ppms[0], ppms[-1]
    w = traces[traces["ppm"] == wide]
    n = traces[traces["ppm"] == narrow]
    m = w.merge(n, on=["file", "target_id", "scan_idx"], suffixes=("_w", "_n"))
    rows = []
    for (file, tid), g in m.groupby(["file", "target_id"], sort=False):
        dw, dn = g["intensity_w"].notna(), g["intensity_n"].notna()
        if not dw.any():
            continue
        both = dw & dn
        rows.append(dict(
            file=file, target_id=tid, ppm_wide=wide, ppm_narrow=narrow,
            scans_detected_wide=int(dw.sum()),
            lost_fraction=float((dw & ~dn).sum() / dw.sum()),
            interference_removed_fraction=float(((g["n_centroids_w"] > 1) & (g["n_centroids_n"] == 1)).mean()),
            peak_changed_fraction=float((both & (g["mz_obs_w"] != g["mz_obs_n"])).sum() / max(both.sum(), 1)),
            mean_intensity_wide=float(g.loc[dw, "intensity_w"].mean()),
            mean_intensity_narrow=float(g.loc[dn, "intensity_n"].mean()) if dn.any() else np.nan,
        ))
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    s = summary[summary["ppm"] == wide][["file", "target_id", "ppm_error_median", "ppm_error_sd"]]
    return out.merge(s, on=["file", "target_id"], how="left")
