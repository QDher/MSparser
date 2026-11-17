"""Data-derived thresholds, window selection, global and per-m/z-group levels."""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.optimize import nnls

MAD = 1.4826
TRACE = ["file", "target_id"]


# ---------------------------------------------------------------- thresholds
def _best_per_trace(win: pd.DataFrame, length: float) -> pd.DataFrame:
    w = win[(win["axis"] == "clock") & win["eligible"] & (win["length"] == length)]
    if w.empty:
        return w
    return w.loc[w.groupby(TRACE)["rsd"].idxmin()]


def rsd_model(mean, thr: dict):
    return np.sqrt(thr["a"] / np.asarray(mean, dtype=float) + thr["b"])


def derive_thresholds(win: pd.DataFrame, cfg: dict) -> dict:
    """RSD: noise model rsd^2 = a/I + b fitted on the most stable window of each
    trace; limit = model x factor. Slope: spread of the slopes of those windows.
    Both are derived at the reference window length and applied to all lengths."""
    tc = cfg["thresholds"]
    lengths = sorted(win.loc[win["eligible"], "length"].unique())
    if not lengths:
        raise ValueError("no eligible window: the acquisitions are too short for the window lengths and minimum scans")
    ref = tc["reference_length_s"] or lengths[len(lengths) // 2]
    ref = min(lengths, key=lambda x: abs(x - ref))
    best = _best_per_trace(win, ref)
    mean, rsd = best["mean"].to_numpy(), np.maximum(best["rsd"].to_numpy(), 1e-6)
    thr = dict(reference_length_s=float(ref), n_traces=int(len(best)), model="noise")
    if len(best) >= 4:
        # relative least squares: minimise sum(((a/I + b) / rsd^2 - 1)^2), a, b >= 0
        A = np.column_stack([1.0 / mean, np.ones_like(mean)]) / (rsd**2)[:, None]
        (a, b), _ = nnls(A, np.ones_like(mean))
    else:
        a, b, thr["model"] = 0.0, float(np.median(rsd**2)), "constant"
    if a == 0 and b == 0:
        b = float(np.median(rsd**2))
    thr.update(a=float(a), b=float(b))
    resid = np.log(rsd / rsd_model(mean, thr))
    spread = MAD * np.median(np.abs(resid - np.median(resid)))
    thr["rsd_factor"] = float(max(np.exp(np.median(resid) + tc["n_sigma"] * spread), tc["min_rsd_factor"]))
    thr["slope_limit"] = float(tc["n_sigma"] * MAD * np.median(np.abs(best["rel_slope"])))
    return thr


def apply_thresholds(win: pd.DataFrame, thr: dict) -> pd.DataFrame:
    win = win.copy()
    win["rsd_limit"] = thr["rsd_factor"] * rsd_model(win["mean"], thr)
    # the slope must be within the limit and estimated precisely enough to be
    # compared with it (short or noisy windows cannot demonstrate stability)
    win["stable"] = (
        win["eligible"] & (win["rsd"] <= win["rsd_limit"])
        & (win["rel_slope"].abs() <= thr["slope_limit"])
        & (win["slope_se"] <= thr["slope_limit"])
    )
    win["sem"] = win["rsd"] / np.sqrt(win["n"].clip(lower=1))
    return win


# ----------------------------------------------------------------- selection
def select_trace_windows(win: pd.DataFrame) -> pd.DataFrame:
    """Per trace and axis: the stable window with the lowest RSD/sqrt(n)."""
    cols = ["start", "length", "n", "mean", "rsd", "rel_slope", "slope_se", "sem"]
    traces = win[["file", "target_id", "axis"]].drop_duplicates()
    st = win[win["stable"]]
    sel = st.loc[st.groupby(TRACE + ["axis"])["sem"].idxmin(), TRACE + ["axis"] + cols] if len(st) else \
        pd.DataFrame(columns=TRACE + ["axis"] + cols)
    out = traces.merge(sel, on=TRACE + ["axis"], how="left")
    out["no_stable_window"] = out["start"].isna()
    return out


def window_grid(win: pd.DataFrame, by: list[str] | None = None) -> pd.DataFrame:
    """Percentage of traces that are stable in each candidate window."""
    by = (by or []) + ["axis"]
    n_all = win.groupby(by)[TRACE].apply(lambda d: len(d.drop_duplicates())).rename("n_traces")
    st = win[win["stable"]]
    n_with = st.groupby(by)[TRACE].apply(lambda d: len(d.drop_duplicates())).rename("n_with_stable")
    g = st.groupby(by + ["length", "start"]).agg(n_stable=("sem", "size"), median_sem=("sem", "median"))
    g = g.reset_index().merge(n_all.reset_index(), on=by).merge(n_with.reset_index(), on=by)
    g["pct_all"] = 100 * g["n_stable"] / g["n_traces"]
    g["pct_with_stable"] = 100 * g["n_stable"] / g["n_with_stable"]
    return g


def best_of_grid(grid: pd.DataFrame, by: list[str] | None = None) -> pd.DataFrame:
    """Highest stable percentage; ties resolved by the lowest median RSD/sqrt(n)."""
    by = (by or []) + ["axis"]
    if grid.empty:
        return grid
    g = grid.sort_values(["pct_all", "median_sem", "start"], ascending=[False, True, True])
    return g.groupby(by, sort=False).head(1).reset_index(drop=True)


def _at_window(win: pd.DataFrame, axis: str, start: float, length: float) -> pd.DataFrame:
    return win[(win["axis"] == axis) & (win["start"] == start) & (win["length"] == length)]


def group_level(win: pd.DataFrame, global_best: pd.DataFrame) -> pd.DataFrame:
    """Per target ion: group window and comparison with the global window."""
    grid = window_grid(win, ["target_id"])
    best = best_of_grid(grid, ["target_id"])
    n_tr = win.groupby(["target_id", "axis"])["file"].nunique().rename("n_acquisitions").reset_index()
    rows = []
    for r in n_tr.itertuples(index=False):
        w = win[(win["target_id"] == r.target_id) & (win["axis"] == r.axis)]
        row = dict(target_id=r.target_id, axis=r.axis, n_acquisitions=r.n_acquisitions)
        b = best[(best["target_id"] == r.target_id) & (best["axis"] == r.axis)]
        gb = global_best[global_best["axis"] == r.axis]
        for name, sel in (("group", b), ("global", gb)):
            if sel.empty:
                row.update({f"{name}_start": np.nan, f"{name}_length": np.nan, f"pct_{name}_window": 0.0})
                continue
            s, length = sel["start"].iloc[0], sel["length"].iloc[0]
            at = _at_window(w, r.axis, s, length)
            el = at[at["eligible"]]
            row.update({
                f"{name}_start": s, f"{name}_length": length,
                f"pct_{name}_window": 100 * at["stable"].sum() / r.n_acquisitions,
                f"median_rsd_{name}_window": el["rsd"].median(),
                f"median_abs_slope_{name}_window": el["rel_slope"].abs().median(),
            })
        row["n_no_stable_window"] = r.n_acquisitions - w.loc[w["stable"], "file"].nunique()
        rows.append(row)
    return pd.DataFrame(rows)


# ------------------------------------------------------------------ profiles
def scaled_profiles(traces: pd.DataFrame, info: pd.DataFrame):
    """Traces scaled by their plateau median on a common time grid (comparison
    only). Returns (long table of scaled traces, summary per time point)."""
    step = float(np.nanmax(info["scan_interval_s"]))
    t_max = traces["time_s"].max()
    grid = np.arange(0, t_max + step / 2, step)
    med = info.set_index(TRACE)["plateau_median"]
    cols = {}
    for key, g in traces.groupby(TRACE, sort=False):
        if key not in med.index or not np.isfinite(med[key]) or med[key] <= 0:
            continue
        g = g.dropna(subset=["intensity"]).sort_values("time_s")
        if len(g) < 2:
            continue
        cols[key] = np.interp(grid, g["time_s"], g["intensity"] / med[key], left=np.nan, right=np.nan)
    if not cols:
        return pd.DataFrame(), pd.DataFrame()
    mat = pd.DataFrame(cols, index=grid)
    summ = pd.DataFrame(dict(
        time_s=grid, n=mat.notna().sum(axis=1).to_numpy(), median=mat.median(axis=1).to_numpy(),
        q25=mat.quantile(0.25, axis=1).to_numpy(), q75=mat.quantile(0.75, axis=1).to_numpy(),
    ))
    long = pd.concat(
        [pd.DataFrame(dict(time_s=grid, file=k[0], target_id=k[1], scaled=v)) for k, v in cols.items()],
        ignore_index=True,
    ).dropna(subset=["scaled"])
    return long, summ


# --------------------------------------------------------------- mixed model
def mixed_models(trace_table: pd.DataFrame) -> str:
    """Dependence of window position and stability on intensity, m/z, reaction
    and run order; acquisition as random effect. Returns a text report."""
    import statsmodels.formula.api as smf

    d = trace_table.copy()
    d["log10_intensity"] = np.log10(d["plateau_median"])
    d["log10_best_rsd"] = np.log10(d["best_rsd"])
    preds = [p for p in ("log10_intensity", "mz", "run_order") if d[p].nunique() > 1]
    if d["reaction"].nunique() > 1:
        preds.append("C(reaction)")
    out = []
    for resp in ("start", "log10_best_rsd", "trend_rel_slope"):
        cols = [resp, "file"] + [p for p in ("log10_intensity", "mz", "run_order", "reaction") if p in d]
        x = d[cols].replace([np.inf, -np.inf], np.nan).dropna()
        out.append(f"=== response: {resp} (n = {len(x)}) ===")
        if not preds or len(x) < len(preds) + 4 or x[resp].nunique() < 2:
            out.append("not fitted: too few traces or no variation\n")
            continue
        formula = f"{resp} ~ " + " + ".join(preds)
        try:
            if x["file"].nunique() < 2 or x["file"].nunique() == len(x):
                raise ValueError("one trace per acquisition: no random effect")
            import warnings
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                fit = smf.mixedlm(formula, x, groups=x["file"]).fit(reml=True)
            if not np.all(np.isfinite(fit.bse)):
                raise ValueError("mixed model did not converge")
            out.append(str(fit.summary()))
        except Exception as e:  # fall back to ordinary least squares
            out.append(f"mixed model not used ({e}); ordinary least squares:")
            try:
                out.append(str(smf.ols(formula, x).fit().summary().tables[1]))
            except Exception as e2:
                out.append(f"not fitted: {e2}")
        out.append("")
    return "\n".join(out)


# ---------------------------------------------------------------- validation
def validate(win_raw: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    """Derive thresholds and windows on a subset of acquisitions; evaluate on the rest."""
    v = cfg["validation"]
    files = np.array(sorted(win_raw["file"].unique()))
    if len(files) < v["min_files"]:
        return pd.DataFrame()
    rng = np.random.default_rng(v["seed"])
    rng.shuffle(files)
    n_test = min(max(1, int(round(len(files) * v["test_fraction"]))), len(files) - 1)
    test_files = set(files[:n_test])
    is_test = win_raw["file"].isin(test_files)
    try:
        thr = derive_thresholds(win_raw[~is_test], cfg)
    except ValueError:
        return pd.DataFrame()
    train, test = apply_thresholds(win_raw[~is_test], thr), apply_thresholds(win_raw[is_test], thr)
    rows = []

    def add(level, best, tr, te):
        for b in best.itertuples(index=False):
            at = _at_window(te, b.axis, b.start, b.length)
            n_te = len(te[te["axis"] == b.axis][TRACE].drop_duplicates())
            rows.append(dict(
                level=level, axis=b.axis, start=b.start, length=b.length,
                n_train=b.n_traces, pct_train=b.pct_all, n_test=n_te,
                pct_test=100 * at["stable"].sum() / n_te if n_te else np.nan,
            ))

    add("global", best_of_grid(window_grid(train)), train, test)
    gbest = best_of_grid(window_grid(train, ["target_id"]), ["target_id"])
    for tid, b in gbest.groupby("target_id"):
        add(tid, b.drop(columns="target_id"), train[train["target_id"] == tid], test[test["target_id"] == tid])
    return pd.DataFrame(rows)


# ------------------------------------------------------- fixed-window fallback
def fixed_window_reproducibility(win: pd.DataFrame, sel: pd.DataFrame, global_best: pd.DataFrame) -> pd.DataFrame:
    """Reproducibility across acquisitions of the mean intensity in the global
    clock-time window, including traces without a stable window."""
    gb = global_best[global_best["axis"] == "clock"]
    if gb.empty:
        return pd.DataFrame()
    at = _at_window(win, "clock", gb["start"].iloc[0], gb["length"].iloc[0])
    at = at.merge(sel[sel["axis"] == "clock"][TRACE + ["no_stable_window"]], on=TRACE)
    rows = []
    for tid, g in at.groupby("target_id"):
        m, ns = g["mean"].dropna(), g.loc[g["no_stable_window"], "mean"].dropna()
        cv = lambda x: float(x.std(ddof=1) / x.mean() * 100) if len(x) > 1 else np.nan
        rows.append(dict(
            target_id=tid, start=gb["start"].iloc[0], length=gb["length"].iloc[0],
            n_acquisitions=len(m), mean_intensity=m.mean() if len(m) else np.nan,
            cv_pct_between_acquisitions=cv(m), n_no_stable_window=len(ns),
            cv_pct_no_stable_window=cv(ns),
        ))
    return pd.DataFrame(rows)
