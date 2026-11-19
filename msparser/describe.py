"""Per-trace description: phases, sliding-window metrics, full-acquisition trend."""
from __future__ import annotations

import numpy as np
import pandas as pd
import ruptures as rpt
from scipy import stats


def segment_phases(t: np.ndarray, y: np.ndarray, cfg: dict) -> dict:
    """Rise / plateau / decay from change points. Undetected scans count as 0."""
    ph = cfg["phases"]
    out = dict(onset_s=np.nan, plateau_start_s=np.nan, plateau_end_s=np.nan, plateau_median=np.nan)
    det = np.isfinite(y)
    if det.sum() < 3:
        return out
    z = np.where(det, y, 0.0)
    n = len(z)
    bkps = [n]
    if n >= 6:
        sigma = 1.4826 * np.median(np.abs(np.diff(z) - np.median(np.diff(z)))) / np.sqrt(2)
        sigma = max(sigma, 1e-9 * max(z.max(), 1.0))
        pen = ph["penalty_factor"] * sigma**2 * np.log(n)
        bkps = rpt.Pelt(model="l2", min_size=3, jump=1).fit(z.reshape(-1, 1)).predict(pen=pen)
    starts = [0] + bkps[:-1]
    means = np.array([z[a:b].mean() for a, b in zip(starts, bkps)])
    k = int(np.argmax(means))
    lo = hi = k
    while lo > 0 and means[lo - 1] >= ph["plateau_fraction"] * means[k]:
        lo -= 1
    while hi < len(means) - 1 and means[hi + 1] >= ph["plateau_fraction"] * means[k]:
        hi += 1
    a, b = starts[lo], bkps[hi]
    vals = y[a:b][det[a:b]]
    if len(vals) == 0:
        return out
    med = float(np.median(vals))
    above = np.flatnonzero(z >= ph["onset_fraction"] * med)
    out.update(
        onset_s=float(t[above[0]]) if len(above) else np.nan,
        plateau_start_s=float(t[a]), plateau_end_s=float(t[b - 1]), plateau_median=med,
    )
    return out


def window_stats(t: np.ndarray, y: np.ndarray, start: float, end: float, w: dict) -> dict:
    """Mean, RSD, relative slope (% per s) and its standard error for the scans
    in [start, end). `eligible`: enough scans, enough of them detected."""
    det = np.isfinite(y)
    a = np.searchsorted(t, start - 1e-9, "left")
    b = np.searchsorted(t, end - 1e-9, "left")
    n_all = b - a
    d = det[a:b]
    n = int(d.sum())
    row = dict(n_scans=n_all, n=n, mean=np.nan, rsd=np.nan, rel_slope=np.nan, slope_se=np.nan)
    if n >= 3:
        tv, yv = t[a:b][d], y[a:b][d]
        mean = yv.mean()
        tc = tv - tv.mean()
        sxx = (tc**2).sum()
        if mean > 0 and sxx > 0:
            slope = (tc * yv).sum() / sxx
            ssr = ((yv - mean - slope * tc) ** 2).sum()
            row.update(mean=mean, rsd=yv.std(ddof=1) / mean, rel_slope=slope / mean * 100.0,
                       slope_se=np.sqrt(ssr / (n - 2) / sxx) / mean * 100.0)
    row["eligible"] = bool(
        n_all >= w["min_scans"] and n_all > 0
        and n / n_all >= w["min_detected_fraction"] and np.isfinite(row["rsd"])
    )
    return row


def sliding_windows(t: np.ndarray, y: np.ndarray, cfg: dict, origin: float = 0.0) -> pd.DataFrame:
    """Metrics for windows [start, start+length) on a grid anchored at `origin`.
    `start` is reported relative to `origin`."""
    w = cfg["windows"]
    step = float(w["step_s"])
    rows = []
    for length in w["lengths_s"]:
        k0 = int(np.ceil((t[0] - origin) / step - 1e-9))
        k1 = int(np.floor((t[-1] - length - origin) / step + 1e-9))
        for k in range(k0, k1 + 1):
            s = origin + k * step
            rows.append(dict(start=round(k * step, 6), length=float(length),
                             **window_stats(t, y, s, s + length, w)))
    return pd.DataFrame(rows)


def full_trend(t: np.ndarray, y: np.ndarray) -> dict:
    """Linear trend over the whole acquisition, in % of mean intensity per second."""
    d = np.isfinite(y)
    out = dict(trend_rel_slope=np.nan, trend_ci_low=np.nan, trend_ci_high=np.nan, trend_p=np.nan)
    if d.sum() < 4 or np.ptp(t[d]) == 0:
        return out
    r = stats.linregress(t[d], y[d])
    mean = y[d].mean()
    half = stats.t.ppf(0.975, d.sum() - 2) * r.stderr
    out.update(
        trend_rel_slope=r.slope / mean * 100, trend_ci_low=(r.slope - half) / mean * 100,
        trend_ci_high=(r.slope + half) / mean * 100, trend_p=r.pvalue,
    )
    return out


def describe_traces(traces: pd.DataFrame, summary: pd.DataFrame, cfg: dict):
    """Return (trace_info, windows). windows has one row per trace, axis, window."""
    keep = summary[summary["analyzed"]][["file", "target_id", "ppm"]]
    info, wins = [], []
    grouped = traces.merge(keep, on=["file", "target_id", "ppm"]).groupby(
        ["file", "target_id", "ppm"], sort=False
    )
    for (file, tid, ppm), g in grouped:
        g = g.sort_values("time_s")
        t, y = g["time_s"].to_numpy(), g["intensity"].to_numpy()
        ph = segment_phases(t, y, cfg)
        info.append(dict(file=file, target_id=tid, ppm=ppm, **ph, **full_trend(t, y),
                         scan_interval_s=float(np.median(np.diff(t))) if len(t) > 1 else np.nan))
        axes = [("clock", 0.0)]
        if np.isfinite(ph["onset_s"]):
            axes.append(("onset", ph["onset_s"]))
        for axis, origin in axes:
            w = sliding_windows(t, y, cfg, origin)
            if len(w):
                w.insert(0, "axis", axis)
                w.insert(0, "ppm", ppm)
                w.insert(0, "target_id", tid)
                w.insert(0, "file", file)
                wins.append(w)
    return pd.DataFrame(info), (pd.concat(wins, ignore_index=True) if wins else pd.DataFrame())
