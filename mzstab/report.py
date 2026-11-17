"""Figures. Tables are written by the CLI."""
from __future__ import annotations

import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap

from .analyze import TRACE, manual_window_for, rsd_model

SURFACE, INK, INK2, MUTED, GRID, AXIS = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
BLUE, ORANGE = "#2a78d6", "#eb6834"
SEQ = LinearSegmentedColormap.from_list("seq_blue", ["#cde2fb", "#86b6ef", "#2a78d6", "#0d366b"])

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "axes.edgecolor": AXIS, "axes.labelcolor": INK2, "axes.titlecolor": INK,
    "axes.titlesize": 11, "axes.titleweight": "bold", "axes.titlelocation": "left",
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True,
    "xtick.color": MUTED, "ytick.color": MUTED, "xtick.labelcolor": INK2, "ytick.labelcolor": INK2,
    "font.size": 9, "lines.linewidth": 2, "legend.frameon": False, "figure.dpi": 150,
})


def _save(fig, path: Path):
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def _shade(ax, start, length, label):
    ax.axvspan(start, start + length, color=BLUE, alpha=0.10, lw=0)
    ax.annotate(label, (start + length / 2, 1.0), xycoords=("data", "axes fraction"),
                ha="center", va="bottom", fontsize=8, color=INK2)


def _outline(ax, start, end, label):
    """Manual window: dashed orange limits, labelled below the automatic label."""
    for v in (start, end):
        ax.axvline(v, color=ORANGE, lw=1.4, ls=(0, (4, 3)))
    ax.annotate(label, ((start + end) / 2, 0.02), xycoords=("data", "axes fraction"),
                ha="center", va="bottom", fontsize=8, color=INK2)


def fig_profiles(long, summ, global_best, path, manual=None):
    fig, ax = plt.subplots(figsize=(7, 4))
    for _, g in long.groupby(TRACE, sort=False):
        ax.plot(g["time_s"], g["scaled"], color=MUTED, lw=0.5, alpha=0.35)
    ax.fill_between(summ["time_s"], summ["q25"], summ["q75"], color=BLUE, alpha=0.25, lw=0,
                    label="Interquartile range")
    ax.plot(summ["time_s"], summ["median"], color=BLUE, label="Median profile")
    gb = global_best[global_best["axis"] == "clock"]
    if len(gb):
        _shade(ax, gb["start"].iloc[0], gb["length"].iloc[0], "global window")
    if manual:
        _outline(ax, float(manual[0]), float(manual[1]), "manual window")
    ax.set(xlabel="Acquisition time (s)", ylabel="Intensity / plateau median",
           ylim=(0, max(1.6, float(np.nanmax(summ["q75"])) * 1.1)))
    ax.set_title("Scaled intensity profiles, all traces", pad=16)
    ax.legend(loc="upper right", ncols=2)
    _save(fig, path)


def fig_heatmap(grid, axis, path):
    g = grid[grid["axis"] == axis]
    if g.empty:
        return
    piv = g.pivot(index="length", columns="start", values="pct_all").sort_index()
    fig, ax = plt.subplots(figsize=(7, 3.6))
    x, y = piv.columns.to_numpy(), np.arange(len(piv))
    dx = np.min(np.diff(x)) if len(x) > 1 else 1.0
    im = ax.pcolormesh(np.append(x, x[-1] + dx), np.append(y, len(y)) - 0.5,
                       np.ma.masked_invalid(piv.to_numpy()), cmap=SEQ, vmin=0, vmax=100,
                       edgecolors=SURFACE, linewidth=0.5)
    ax.set_yticks(y, [f"{v:g}" for v in piv.index])
    ax.grid(False)
    ax.set(xlabel=("Window start, acquisition time (s)" if axis == "clock"
                   else "Window start, time after signal onset (s)"),
           ylabel="Window length (s)", title="Traces with a stable signal in each window")
    fig.colorbar(im, ax=ax, label="% of traces")
    _save(fig, path)


def fig_rsd_intensity(win, thr, path):
    w = win[(win["axis"] == "clock") & win["eligible"] & (win["length"] == thr["reference_length_s"])]
    if w.empty:
        return
    best = w.loc[w.groupby(TRACE)["rsd"].idxmin()]
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.scatter(best["mean"], best["rsd"] * 100, s=18, color=MUTED, alpha=0.7, lw=0,
               label="Most stable window of each trace")
    xs = np.geomspace(best["mean"].min(), best["mean"].max(), 100)
    ax.plot(xs, rsd_model(xs, thr) * 100, color=BLUE, label="Noise model")
    ax.plot(xs, rsd_model(xs, thr) * 100 * thr["rsd_factor"], color=ORANGE, label="Stability limit")
    ax.set(xscale="log", yscale="log", xlabel="Mean intensity in window", ylabel="RSD (%)",
           title=f"RSD versus intensity ({thr['reference_length_s']:g} s windows)")
    ax.legend()
    _save(fig, path)


def fig_threshold_hists(win, thr, path):
    w = win[(win["axis"] == "clock") & win["eligible"] & (win["length"] == thr["reference_length_s"])]
    if w.empty:
        return
    lim = thr["slope_limit"]
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.6))
    ratio = w["rsd"] / rsd_model(w["mean"], thr)
    axes[0].hist(np.log10(ratio), bins=40, color=BLUE, edgecolor=SURFACE, lw=0.5)
    axes[0].axvline(np.log10(thr["rsd_factor"]), color=ORANGE)
    axes[0].set(xlabel="log10(RSD / noise model)", ylabel="Windows", title="RSD relative to noise model")
    s = w["rel_slope"].clip(*np.nanpercentile(w["rel_slope"], [1, 99]))
    axes[1].hist(s, bins=40, color=BLUE, edgecolor=SURFACE, lw=0.5)
    for v in (-lim, lim):
        axes[1].axvline(v, color=ORANGE)
    axes[1].set(xlabel="Relative slope (% per s)", ylabel="Windows", title="Slope (limits in orange)")
    _save(fig, path)


def fig_mass_error(traces, ppms, path):
    wide = max(ppms)
    e = traces.loc[traces["ppm"] == wide, "ppm_error"].dropna()
    if e.empty:
        return
    fig, ax = plt.subplots(figsize=(6, 3.6))
    ax.hist(e, bins=60, range=(-wide, wide), color=BLUE, edgecolor=SURFACE, lw=0.5)
    for p in ppms:
        if p < wide:
            for v in (-p, p):
                ax.axvline(v, color=ORANGE)
            ax.annotate(f"±{p:g} ppm", (p, 0.95), xycoords=("data", "axes fraction"),
                        ha="left", va="top", fontsize=8, color=INK2, xytext=(3, 0),
                        textcoords="offset points")
    ax.set(xlabel="Mass error (ppm)", ylabel="Scans",
           title=f"Mass error of detected targets (±{wide:g} ppm extraction)")
    _save(fig, path)


def fig_target_traces(traces, sel, groups, out_dir: Path, manual=None):
    out_dir.mkdir(parents=True, exist_ok=True)
    sel = sel[sel["axis"] == "clock"].set_index(TRACE)
    for tid, g in traces.groupby("target_id", sort=False):
        fig, ax = plt.subplots(figsize=(7, 3.6))
        for file, t in g.groupby("file", sort=False):
            t = t.sort_values("time_s")
            ax.plot(t["time_s"], t["intensity"], color=MUTED, lw=0.8, alpha=0.7)
            if (file, tid) in sel.index and not sel.loc[(file, tid), "no_stable_window"]:
                s = sel.loc[(file, tid)]
                m = (t["time_s"] >= s["start"]) & (t["time_s"] < s["start"] + s["length"])
                ax.plot(t.loc[m, "time_s"], t.loc[m, "intensity"], color=BLUE, lw=1.6)
        gr = groups[(groups["target_id"] == tid) & (groups["axis"] == "clock")]
        if len(gr) and np.isfinite(gr["group_start"].iloc[0]):
            _shade(ax, gr["group_start"].iloc[0], gr["group_length"].iloc[0], "group window")
        mw = manual_window_for(tid, manual)
        if mw:
            _outline(ax, mw[0], mw[1], "manual window")
        ax.set(xlabel="Acquisition time (s)", ylabel="Absolute intensity", ylim=(0, None))
        ax.set_title(f"{tid}: traces per acquisition (blue = selected stable window)", pad=16)
        _save(fig, out_dir / (re.sub(r"[^\w.+-]+", "_", tid) + ".png"))
