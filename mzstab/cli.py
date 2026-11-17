"""Command line: convert -> targets -> ingest -> extract -> analyze (+ report)."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from . import analyze as an
from . import report as rp
from .config import load_config
from .convert import convert_all
from .describe import describe_traces
from .detect import load_detected, resolve_resolution, resolve_windows, save_detected
from .extract import compare_windows, extract_file
from .io import find_input, ingest_file, load_samples
from .targets import generate_targets, load_components


def _work(cfg) -> Path:
    p = Path(cfg["paths"]["work_dir"])
    p.mkdir(parents=True, exist_ok=True)
    return p


def cmd_convert(cfg, args):
    cmds = convert_all(cfg, dry_run=args.dry_run)
    if args.dry_run:
        for c in cmds:
            print(" ".join(c))


def cmd_targets(cfg, args=None):
    """Write the proposed target table; edit the 'use' column to confirm."""
    out = _work(cfg) / "targets.csv"
    if out.exists() and not getattr(args, "force", False):
        print(f"{out} exists (kept; use --force to regenerate)")
        return
    t = generate_targets(load_components(cfg["paths"]["components"]), cfg["targets"]["pairs"])
    t.to_csv(out, index=False)
    print(f"{len(t)} targets written to {out}: {t['type'].value_counts().to_dict()}")


def cmd_ingest(cfg, args=None):
    work, mzml_dir = _work(cfg), Path(cfg["paths"]["mzml_dir"])
    dropped, detected = [], {}
    for f in load_samples(cfg)["file"]:
        path = find_input(mzml_dir, f)
        s, detected[f] = ingest_file(path, cfg, work)
        d = s[~s["kept"]].groupby("drop_reason").size().rename("n_scans").reset_index()
        d.insert(0, "file", f)
        dropped.append(d)
        d = detected[f]
        print(f"{f}: {int(s['kept'].sum())}/{len(s)} scans kept | {d['vendor']}, "
              f"{'+' if d['polarity'] > 0 else '-'}, {d['duration_s']:.0f} s, "
              f"{1 / d['scan_interval_s']:.1f} scans/s" if d["n_scans"] > 1 else f"{f}: no usable scan")
    save_detected(work, detected)
    pd.concat(dropped, ignore_index=True).to_csv(work / "scans_dropped.csv", index=False)


def cmd_extract(cfg, args=None):
    work = _work(cfg)
    if not (work / "targets.csv").exists():
        cmd_targets(cfg)
    targets = pd.read_csv(work / "targets.csv")
    known = set(targets["component"])
    traces, summary = [], []
    detected = load_detected(work)
    for s in load_samples(cfg).itertuples(index=False):
        injected = {c.strip() for c in s.components.split(";") if c.strip()} or None
        if injected and not injected <= known:
            raise ValueError(f"{s.file}: unknown components {sorted(injected - known)}")
        t, m = extract_file(s.file, cfg, targets, injected, work,
                            resolve_resolution(cfg, detected.get(s.file)))
        traces.append(t)
        summary.append(m)
    traces = pd.concat(traces, ignore_index=True)
    summary = pd.concat(summary, ignore_index=True)
    traces.to_parquet(work / "traces.parquet", index=False)
    summary.to_parquet(work / "trace_summary.parquet", index=False)
    compare_windows(traces, summary).to_parquet(work / "window_comparison.parquet", index=False)
    n = summary.drop_duplicates(["file", "target_id"])
    print(f"{len(n)} traces; analyzed per ppm setting: "
          f"{summary.groupby('ppm')['analyzed'].sum().astype(int).to_dict()}")


def cmd_analyze(cfg, args=None):
    work, rep = _work(cfg), Path(cfg["paths"]["report_dir"])
    rep.mkdir(parents=True, exist_ok=True)
    traces = pd.read_parquet(work / "traces.parquet")
    summary = pd.read_parquet(work / "trace_summary.parquet")
    comparison = pd.read_parquet(work / "window_comparison.parquet")
    samples = load_samples(cfg)[["file", "reaction", "run_order"]]
    detected = load_detected(work)
    used = [f for f in detected if f in set(samples["file"])]
    cfg["windows"] = resolve_windows(cfg, {f: detected[f] for f in used})
    w = cfg["windows"]
    print(f"time windows: lengths {w['lengths_s']} s, step {w['step_s']:g} s, "
          f"at least {w['min_scans']} scans" + (" (automatic)" if w["source"] else ""))
    (rep / "settings_used.json").write_text(json.dumps({
        "ppm_windows": cfg["extraction"]["ppm_windows"], "windows": w,
        "files": {f: {**detected[f], "resolution": resolve_resolution(cfg, detected[f])} for f in used},
    }, indent=2))
    info, windows = describe_traces(traces, summary, cfg)
    if windows.empty:
        raise SystemExit("no analyzable trace: see trace_summary for flags")
    ppms = cfg["extraction"]["ppm_windows"]
    tables = {k: [] for k in ("trace_results", "global_grid", "global_windows", "group_windows",
                              "validation", "fixed_window", "window_choice", "manual_window_traces",
                              "file_summary", "threshold_points")}
    thresholds = {}
    for ppm in ppms:
        tag = f"ppm{ppm:g}"
        w_raw = windows[windows["ppm"] == ppm].drop(columns="ppm")
        inf = info[info["ppm"] == ppm].drop(columns="ppm")
        tr = traces[traces["ppm"] == ppm].merge(inf[an.TRACE], on=an.TRACE)
        if w_raw.empty:
            continue
        thr = an.derive_thresholds(w_raw, cfg)
        an.set_limits(thr, cfg)
        if thr["n_traces"] < 5:
            print(f"ppm {ppm:g}: the noise model rests on only {thr['n_traces']} trace(s) and is not reliable. "
                  "Check trace_summary for undetected targets, or add acquisitions.")
        print(f"ppm {ppm:g}: RSD limit = {thr['rsd_factor']:g} x noise model ({thr['rsd_factor_source']}); "
              f"slope limit = {thr['slope_limit']:g} % per s ({thr['slope_limit_source']})")
        thresholds[tag] = thr
        win = an.apply_thresholds(w_raw, thr)
        sel = an.select_trace_windows(win)
        grid = an.window_grid(win)
        gbest = an.best_of_grid(grid)
        groups = an.group_level(win, gbest)

        el = win[(win["axis"] == "clock") & win["eligible"]]
        best_rsd = el.groupby(an.TRACE)["rsd"].min().rename("best_rsd").reset_index()
        res = (sel.merge(inf, on=an.TRACE).merge(best_rsd, on=an.TRACE, how="left")
               .merge(summary[summary["ppm"] == ppm].drop(columns="ppm"), on=an.TRACE)
               .merge(samples, on="file", how="left"))
        (rep / f"mixed_model_{tag}.txt").write_text(an.mixed_models(res[res["axis"] == "clock"]))

        manual = cfg.get("manual_windows") or {}
        choice, manual_tr = an.window_choice(tr, thr, gbest, groups, manual, cfg)
        # the points behind the RSD-versus-intensity figure: best window of each trace
        pts = an._best_per_trace(w_raw, thr["reference_length_s"])[
            an.TRACE + ["start", "length", "n", "mean", "rsd", "rel_slope"]]
        pts = pts.merge(summary.loc[summary["ppm"] == ppm, an.TRACE + ["mz"]], on=an.TRACE, how="left")
        per_file = an.file_summary(win, sel, gbest).merge(samples[["file", "reaction"]], on="file", how="left")
        for name, df in (("window_choice", choice), ("manual_window_traces", manual_tr), ("file_summary", per_file), ("threshold_points", pts),
                         ("trace_results", res), ("global_grid", grid), ("global_windows", gbest),
                         ("group_windows", groups), ("validation", an.validate(w_raw, cfg)),
                         ("fixed_window", an.fixed_window_reproducibility(win, sel, gbest))):
            if len(df):
                tables[name].append(df.assign(ppm=ppm))

        long, summ = an.scaled_profiles(tr, inf)
        if len(summ):
            summ.to_csv(rep / f"profile_{tag}.csv", index=False)
            rp.fig_profiles(long, summ, gbest, rep / f"profiles_{tag}.png", manual.get("global"))
        for axis in ("clock", "onset"):
            rp.fig_heatmap(grid, axis, rep / f"stable_pct_{axis}_{tag}.png")
        rp.fig_rsd_intensity(win, thr, rep / f"rsd_vs_intensity_{tag}.png")
        rp.fig_threshold_hists(win, thr, rep / f"threshold_histograms_{tag}.png")
        rp.fig_target_traces(tr, sel, groups, rep / f"traces_{tag}", manual)

    for name, parts in tables.items():
        if parts:
            pd.concat(parts, ignore_index=True).to_csv(rep / f"{name}.csv", index=False)
        elif (rep / f"{name}.csv").exists():
            (rep / f"{name}.csv").unlink()  # e.g. manual windows removed since the last run
    (rep / "thresholds.json").write_text(json.dumps(thresholds, indent=2))
    summary.to_csv(rep / "trace_summary.csv", index=False)
    pd.read_csv(work / "targets.csv").to_csv(rep / "targets.csv", index=False)
    if (work / "scans_dropped.csv").exists():
        pd.read_csv(work / "scans_dropped.csv").to_csv(rep / "scans_dropped.csv", index=False)
    rp.fig_mass_error(traces, ppms, rep / "mass_error.png")

    # effect of the m/z window on the stability results
    if len(comparison) and tables["trace_results"]:
        res = pd.concat(tables["trace_results"], ignore_index=True)
        res = res[res["axis"] == "clock"]
        for label, ppm in (("wide", max(ppms)), ("narrow", min(ppms))):
            r = res[res["ppm"] == ppm][an.TRACE + ["best_rsd", "rsd", "rel_slope", "no_stable_window"]]
            r = r.rename(columns={"best_rsd": f"best_rsd_{label}", "rsd": f"selected_rsd_{label}",
                                  "rel_slope": f"selected_slope_{label}",
                                  "no_stable_window": f"no_stable_window_{label}"})
            comparison = comparison.merge(r, on=an.TRACE, how="left")
        comparison.to_csv(rep / "window_comparison.csv", index=False)

    gw = pd.concat(tables["global_windows"]) if tables["global_windows"] else pd.DataFrame()
    if gw.empty:
        print("no window is stable in any trace: no window is recommended")
    for r in gw.itertuples(index=False):
        ref = "acquisition time" if r.axis == "clock" else "time after signal onset"
        print(f"ppm {r.ppm:g} | {ref}: global window {r.start:g}-{r.start + r.length:g} s, "
              f"stable in {r.pct_all:.0f}% of {r.n_traces} traces")
    print(f"report written to {rep}")


def cmd_run(cfg, args=None):
    cmd_targets(cfg)
    cmd_ingest(cfg)
    cmd_extract(cfg)
    cmd_analyze(cfg)


def cmd_ui(args):
    import subprocess
    import sys

    app = Path(__file__).with_name("app.py")
    subprocess.run([sys.executable, "-m", "streamlit", "run", str(app), "--", args.project])


def main(argv=None):
    p = argparse.ArgumentParser(prog="mzstab", description=__doc__)
    sub = p.add_subparsers(dest="cmd", required=True)
    for name, fn, hlp in (
        ("convert", cmd_convert, "vendor raw -> centroided mzML"),
        ("targets", cmd_targets, "propose target ions from components"),
        ("ingest", cmd_ingest, "mzML -> Parquet, scan filtering"),
        ("extract", cmd_extract, "extract ion traces and QC flags"),
        ("analyze", cmd_analyze, "thresholds, window selection, report"),
        ("run", cmd_run, "targets + ingest + extract + analyze"),
    ):
        s = sub.add_parser(name, help=hlp)
        s.add_argument("-c", "--config", required=True)
        s.set_defaults(fn=fn)
        if name == "convert":
            s.add_argument("--dry-run", action="store_true", help="print commands only")
        if name == "targets":
            s.add_argument("--force", action="store_true", help="overwrite targets.csv")
    u = sub.add_parser("ui", help="graphical interface in the browser (needs streamlit)")
    u.add_argument("project", nargs="?", default=".", help="project folder")
    args = p.parse_args(argv)
    if args.cmd == "ui":
        return cmd_ui(args)
    args.fn(load_config(args.config), args)


if __name__ == "__main__":
    main()
