"""Streamlit interface. Start with:  mzstab ui <project folder>

A project folder holds one vendor run: config.yaml, components.yaml,
samples.csv, mzml/, work/, report/. The interface writes the same files the
command line uses, so every run can be repeated with `mzstab run -c config.yaml`.
"""
from __future__ import annotations

import argparse
import contextlib
import copy
import io
import json
import shutil
import sys
from pathlib import Path

import altair as alt
import numpy as np
import pandas as pd
import streamlit as st
import yaml

from mzstab import cli
from mzstab.config import DEFAULTS, _merge, load_config
from mzstab.detect import load_detected, resolve_resolution, resolve_windows
from mzstab.analyze import manual_window_for
from mzstab.targets import Component, generate_targets, normalize_mode

CHARGE_COLS = {"+1": 1, "+2": 2, "-1": -1, "-2": -2}
# labels shown in the interface for the two charge modes stored in components.yaml
MODE_LABELS = {"proton": "Molecule (gains or loses H⁺)", "intrinsic": "Ion (already charged)"}
MODE_VALUES = {v: k for k, v in MODE_LABELS.items()}
GREY, BLUE = "#898781", "#2a78d6"

st.set_page_config(page_title="mzstab", layout="wide")


# ------------------------------------------------------------------ project
def user_cfg(root: Path) -> dict:
    p = root / "config.yaml"
    user = (yaml.safe_load(p.read_text()) or {}) if p.exists() else {}
    return _merge(copy.deepcopy(DEFAULTS), user)


def save_cfg(root: Path, cfg: dict):
    (root / "config.yaml").write_text(yaml.safe_dump(cfg, sort_keys=False))


def resolved(root: Path) -> dict:
    return load_config(root / "config.yaml")


def parse_numbers(text: str) -> list[float]:
    vals = [float(x) for x in text.replace(";", ",").split(",") if x.strip()]
    if not vals or min(vals) <= 0:
        raise ValueError("enter positive numbers separated by commas")
    return vals


def opt_float(text: str):
    return float(text) if str(text).strip() else None


def run_steps(root: Path, steps):
    """Run pipeline steps, capturing printed output into the session log."""
    buf = io.StringIO()
    ok = True
    for label, fn, args in steps:
        try:
            with st.spinner(f"{label}..."), contextlib.redirect_stdout(buf):
                fn(resolved(root), args)
            buf.write(f"[{label}: done]\n")
        except (Exception, SystemExit) as e:
            buf.write(f"[{label}: FAILED] {e}\n")
            st.error(f"{label} failed: {e}")
            ok = False
            break
    st.session_state["log"] = buf.getvalue()
    st.session_state["refresh"] = ok
    return ok


# ------------------------------------------------------------------ sidebar
default_root = sys.argv[1] if len(sys.argv) > 1 else "."
st.sidebar.title("mzstab")
root = Path(st.sidebar.text_input("Project folder", str(Path(default_root).resolve()))).expanduser()
if not root.exists():
    st.sidebar.warning("Folder does not exist.")
    if st.sidebar.button("Create folder"):
        root.mkdir(parents=True)
        st.rerun()
    st.stop()
if not (root / "config.yaml").exists():
    save_cfg(root, copy.deepcopy(DEFAULTS))
cfg = user_cfg(root)
paths = resolved(root)["paths"]
work, rep = Path(paths["work_dir"]), Path(paths["report_dir"])
status = {
    "Components": Path(paths["components"]).exists(),
    "Samples": Path(paths["samples"]).exists(),
    "Targets": (work / "targets.csv").exists(),
    "Traces": (work / "traces.parquet").exists(),
    "Report": (rep / "trace_results.csv").exists(),
}
for name, done in status.items():
    st.sidebar.write(("✅ " if done else "⬜ ") + name)

tabs = st.tabs(["1. Settings", "2. Components", "3. Samples", "4. Targets", "5. Run", "6. Results"])

# ----------------------------------------------------------------- settings
def _txt(v) -> str:
    return "" if v in (None, "auto") else (", ".join(f"{x:g}" for x in v) if isinstance(v, list) else str(v))


with tabs[0]:
    st.caption("Only the m/z windows need to be chosen. The other settings are read from the data files "
               "or use default values; they can be changed under 'Advanced settings'.")
    detected = load_detected(work)
    with st.form("settings"):
        ppm = st.text_input("m/z windows, ± ppm (all are run and compared)",
                            ", ".join(f"{p:g}" for p in cfg["extraction"]["ppm_windows"]))
        with st.expander("Advanced settings (leave empty for automatic)"):
            c1, c2, c3 = st.columns(3)
            with c1:
                st.markdown("**Instrument**")
                vendors = ["auto", "thermo", "waters"]
                vendor = st.selectbox("Vendor", vendors, index=vendors.index(cfg["vendor"]),
                                      help="auto: recognised from the files")
                res_value = st.text_input("Resolving power", _txt(cfg["resolution"].get("value")),
                                          help="Empty: read from the file if present, otherwise "
                                               "60,000 at m/z 200 (Orbitrap) or 20,000 (TOF). "
                                               "Only used to flag targets too close to be separated.")
                res_at = st.number_input("Specified at m/z (Orbitrap)", value=float(cfg["resolution"].get("at_mz", 200)))
                max_it = st.text_input("Maximum injection time, ms (Orbitrap)", _txt(cfg["extraction"]["max_injection_time_ms"]),
                                       help="Empty: the largest injection time observed in the file is used")
                sat = st.text_input("Saturation intensity (TOF)", _txt(cfg["extraction"]["saturation_intensity"]),
                                    help="Empty: the saturation flag is not evaluated")
                functions = st.text_input("Waters functions to keep", _txt(cfg["scan_filter"]["functions"]),
                                          help="Empty: function 1")
            with c2:
                st.markdown("**Time windows**")
                lengths = st.text_input("Window lengths, s", _txt(cfg["windows"]["lengths_s"]),
                                        help="Empty: chosen from the scan rate and the acquisition duration")
                step_s = st.text_input("Step between window starts, s", _txt(cfg["windows"]["step_s"]),
                                       help="Empty: 1 s, or the scan interval if longer")
                min_scans = st.number_input("Minimum scans per window", value=int(cfg["windows"]["min_scans"]),
                                            min_value=3, step=1)
                st.markdown("**m/z extraction**")
                min_mda = st.number_input("Minimum half-width, mDa", value=float(cfg["extraction"]["min_window_mda"]),
                                          step=0.1, format="%.2f")
                min_run = st.number_input("Detection: consecutive scans", value=int(cfg["extraction"]["min_consecutive_scans"]),
                                          min_value=1, step=1)
            with c3:
                st.markdown("**Folders**")
                mzml_dir = st.text_input("Data folder (mzML or CSV files)", cfg["paths"]["mzml_dir"])
                raw_dir = st.text_input("Raw folder (for conversion)", cfg["paths"]["raw_dir"])
                st.markdown("**Converters**")
                thermo_parser = st.text_input("ThermoRawFileParser program or its folder", cfg["convert"]["thermo_parser"],
                                              help="Program name if it is on the PATH, otherwise its full path, "
                                                   "e.g. C:\\Tools\\ThermoRawFileParser\\ThermoRawFileParser.exe")
                msconvert = st.text_input("msconvert program (Waters)", cfg["convert"]["msconvert"],
                                          help="Full path of msconvert.exe if it is not on the PATH")
                st.markdown("**Conversion (Waters lock mass)**")
                lock = st.text_input("Lock-mass m/z, positive", _txt(cfg["convert"]["lockmass_mz"]))
                lock_neg = st.text_input("Lock-mass m/z, negative", _txt(cfg["convert"]["lockmass_mz_neg"]))
                st.markdown("**CSV input**")
                csv_unit = st.selectbox("Time unit in CSV files", ["s", "min"],
                                        index=["s", "min"].index(cfg["csv"]["time_unit"]))
                csv_pol = st.selectbox("Polarity of CSV files", [1, -1], index=[1, -1].index(int(cfg["csv"]["polarity"])),
                                       format_func={1: "positive", -1: "negative"}.get)
        if st.form_submit_button("Save settings", type="primary"):
            try:
                cfg["vendor"] = vendor
                cfg["resolution"] = {"value": opt_float(res_value) or "auto", "at_mz": res_at}
                cfg["extraction"].update(
                    ppm_windows=sorted(parse_numbers(ppm), reverse=True), min_window_mda=min_mda,
                    min_consecutive_scans=int(min_run), max_injection_time_ms=opt_float(max_it),
                    saturation_intensity=opt_float(sat))
                cfg["scan_filter"]["functions"] = [int(x) for x in parse_numbers(functions)] if functions.strip() else None
                cfg["paths"].update(mzml_dir=mzml_dir, raw_dir=raw_dir)
                cfg["windows"].update(lengths_s=sorted(parse_numbers(lengths)) if lengths.strip() else "auto",
                                      step_s=opt_float(step_s) or "auto", min_scans=int(min_scans))
                cfg["convert"].update(lockmass_mz=opt_float(lock), lockmass_mz_neg=opt_float(lock_neg),
                                      thermo_parser=thermo_parser.strip() or "ThermoRawFileParser",
                                      msconvert=msconvert.strip() or "msconvert")
                cfg["csv"].update(time_unit=csv_unit, polarity=int(csv_pol))
                save_cfg(root, cfg)
                st.success("Settings saved.")
            except ValueError as e:
                st.error(f"Not saved: {e}")

    st.subheader("Read from the data files")
    if not detected:
        st.info("Shown after the files have been read (tab 5, 'Run analysis').")
    else:
        rcfg = resolved(root)
        rows = []
        for f, d in detected.items():
            r = resolve_resolution(rcfg, d)
            rows.append({
                "file": f, "vendor": d["vendor"], "polarity": "+" if d["polarity"] > 0 else "-",
                "scan range": "–".join(f"{x:g}" for x in d["scan_range"]) if np.isfinite(d["scan_range"][0]) else "not in file",
                "scans": d["n_scans"], "dropped scans": d["n_dropped"],
                "duration (s)": round(d["duration_s"], 1),
                "scans per s": round(1 / d["scan_interval_s"], 2) if d["n_scans"] > 1 else None,
                "resolving power": f"{r['value']:,.0f} ({r['source']})",
            })
        st.dataframe(pd.DataFrame(rows), hide_index=True, width="stretch")
        w = resolve_windows(rcfg, detected)
        st.caption(f"Time windows used: lengths {', '.join(f'{x:g}' for x in w['lengths_s'])} s; "
                   f"step {w['step_s']:g} s; at least {w['min_scans']} scans per window.")

# --------------------------------------------------------------- components
with tabs[1]:
    st.caption("Enter every injected molecule or ion as a neutral formula (for example Na, not Na⁺), tick the "
               "charge states that are possible, and choose what the component is:")
    st.markdown(
        f"- **{MODE_LABELS['proton']}**: a neutral molecule that is protonated or deprotonated in the source, "
        "e.g. [M+H]⁺, [M−H]⁻.\n"
        f"- **{MODE_LABELS['intrinsic']}**: a species that is an ion by itself (Na, K, NH4, Cl, a quaternary "
        "ammonium…). No H⁺ is added or removed. **Choose this for an ion to obtain its adducts with the "
        "molecules**: sodium with +1 gives [M+Na]⁺; chloride with −1 gives [M+Cl]⁻.\n"
        "- Ions made of two different components are also proposed: two molecules sharing H⁺ "
        "([M1+M2+H]⁺), or a cation with an anion.")
    cp = Path(paths["components"])
    rows = []
    if cp.exists():
        for d in (yaml.safe_load(cp.read_text()) or {}).get("components", []):
            rows.append({"name": d["name"], "formula": d["formula"],
                         **{k: v in d.get("charges", []) for k, v in CHARGE_COLS.items()},
                         "mode": MODE_LABELS[normalize_mode(d.get("mode", "proton"))]})
    comp_df = pd.DataFrame(rows, columns=["name", "formula", *CHARGE_COLS, "mode"])
    for k in CHARGE_COLS:
        comp_df[k] = comp_df[k].astype(bool)
    edited = st.data_editor(
        comp_df, num_rows="dynamic", width="stretch", key="components_editor",
        column_config={
            "name": st.column_config.TextColumn("Name", required=True),
            "formula": st.column_config.TextColumn("Neutral formula", required=True),
            **{k: st.column_config.CheckboxColumn(k, default=False) for k in CHARGE_COLS},
            "mode": st.column_config.SelectboxColumn("The component is a…", options=list(MODE_VALUES),
                                                     default=MODE_LABELS["proton"], required=True, width="medium"),
        })
    if st.button("Save components", type="primary"):
        try:
            out, comps = [], []
            for r in edited.dropna(subset=["name", "formula"]).to_dict("records"):
                charges = [v for k, v in CHARGE_COLS.items() if r.get(k) is True]
                name, mode = str(r["name"]).strip(), MODE_VALUES.get(r.get("mode"), "proton")
                if not charges:
                    raise ValueError(f"{name}: tick at least one charge state")
                c = Component(name, str(r["formula"]).strip(), tuple(charges), mode)
                c.composition()  # validates the formula
                comps.append(c)
                out.append({"name": name, "formula": c.formula, "charges": charges, "mode": mode})
            if len({c.name for c in comps}) != len(comps):
                raise ValueError("component names must be unique")
            generate_targets(comps, cfg["targets"]["pairs"])
            cp.write_text(yaml.safe_dump({"components": out}, sort_keys=False))
            st.success(f"{len(out)} components saved. Regenerate the targets in tab 4.")
        except Exception as e:
            st.error(f"Not saved: {e}")

# ------------------------------------------------------------------ samples
with tabs[2]:
    st.caption("One row per acquisition. 'file' is the file name without extension. 'components' lists what "
               "was injected, separated by ';' (empty = all components).")
    sp = Path(paths["samples"])
    cols = ["file", "reaction", "components", "run_order"]
    if "samples_df" not in st.session_state or st.session_state.get("samples_root") != str(root):
        st.session_state["samples_df"] = (
            pd.read_csv(sp, dtype=str, keep_default_na=False) if sp.exists() else pd.DataFrame(columns=cols)
        ).reindex(columns=cols).fillna("")
        st.session_state["samples_root"] = str(root)
    if st.button("Add files found in the data and raw folders"):
        found = sorted({p.stem for d in (paths["mzml_dir"], paths["raw_dir"]) if Path(d).exists()
                        for p in Path(d).iterdir() if p.suffix.lower() in (".mzml", ".csv", ".raw")})
        cur = st.session_state["samples_df"]
        new = [f for f in found if f not in set(cur["file"])]
        add = pd.DataFrame({"file": new, "reaction": "", "components": "",
                            "run_order": [str(len(cur) + i + 1) for i in range(len(new))]})
        st.session_state["samples_df"] = pd.concat([cur, add], ignore_index=True)
        st.session_state.pop("samples_editor", None)
        st.info(f"{len(new)} file(s) added." if new else "No new file found.")
    s_edit = st.data_editor(st.session_state["samples_df"], num_rows="dynamic", width="stretch",
                            key="samples_editor",
                            column_config={"file": st.column_config.TextColumn("file", required=True)})
    if st.button("Save samples", type="primary"):
        s_out = s_edit.fillna("").astype(str)
        s_out = s_out[s_out["file"].str.strip() != ""]
        if s_out["file"].duplicated().any():
            st.error("Not saved: duplicate file names.")
        else:
            s_out.to_csv(sp, index=False)
            st.session_state["samples_df"] = s_out.reset_index(drop=True)
            st.success(f"{len(s_out)} acquisitions saved.")

# ------------------------------------------------------------------ targets
with tabs[3]:
    st.caption("Ions proposed from the components. Untick 'use' for ions that should not be followed.")
    tp = work / "targets.csv"
    if st.button("Generate targets from components" if not tp.exists()
                 else "Regenerate targets (discards the current selection)"):
        run_steps(root, [("Targets", cli.cmd_targets, argparse.Namespace(force=True))])
    if tp.exists():
        t = pd.read_csv(tp)
        t["use"] = t["use"].astype(bool)
        if "type" in t:
            n = t["type"].value_counts()
            st.caption(f"{len(t)} ions: {n.get('ion', 0)} of single components, {n.get('adduct', 0)} adducts of a "
                       f"molecule with an entered ion, {n.get('pair', 0)} made of two components. Positive ions "
                       "are followed in positive-mode acquisitions, negative ions in negative-mode ones.")
            cps = (yaml.safe_load(Path(paths["components"]).read_text()) or {}).get("components", []) \
                if Path(paths["components"]).exists() else []
            if cps and not any(normalize_mode(c.get("mode", "proton")) == "intrinsic" for c in cps):
                st.warning(f"No component is set as '{MODE_LABELS['intrinsic']}', so no adduct with an entered ion "
                           "(such as [M+Na]⁺) is proposed. If a component is an ion (Na, K, Cl…), set its "
                           "type to this in tab 2 and regenerate the targets.")
        t_edit = st.data_editor(
            t, width="stretch", hide_index=True, key="targets_editor",
            disabled=[c for c in t.columns if c != "use"],
            column_order=["use", "type", "component", "ion", "charge", "mz", "iso_mz", "iso_ratio", "requires"],
            column_config={"mz": st.column_config.NumberColumn("m/z", format="%.5f"),
                           "iso_mz": st.column_config.NumberColumn("isotope m/z", format="%.5f"),
                           "iso_ratio": st.column_config.NumberColumn("isotope ratio", format="%.4f"),
                           "requires": "requires components"})
        if st.button("Save target selection", type="primary"):
            t_edit.assign(use=t_edit["use"].astype(int)).to_csv(tp, index=False)
            st.success(f"{int(t_edit['use'].sum())} of {len(t_edit)} targets selected.")

# ---------------------------------------------------------------------- run
with tabs[4]:
    ready = status["Components"] and status["Samples"]
    if not ready:
        st.warning("Save the components (tab 2) and the samples (tab 3) first.")
    st.caption("Conversion needs ThermoRawFileParser or msconvert installed (program path under Settings > "
               "Advanced settings > Converters). It is skipped for files that already have an mzML. If mzML "
               "files are already in the data folder, start at 'Run analysis'.")
    st.caption(f"Raw folder: {paths['raw_dir']} — data folder: {paths['mzml_dir']}")
    b = st.columns(5)
    none = argparse.Namespace(force=False, dry_run=False)
    if b[0].button("Convert raw files", disabled=not ready):
        run_steps(root, [("Convert", cli.cmd_convert, none)])
    if b[1].button("Run analysis", type="primary", disabled=not ready,
                   help="Read mzML, extract traces, derive thresholds, select windows"):
        run_steps(root, [("Targets", cli.cmd_targets, none), ("Read mzML", cli.cmd_ingest, none),
                         ("Extract traces", cli.cmd_extract, none), ("Analyze", cli.cmd_analyze, none)])
    if b[2].button("Re-extract only", disabled=not (work / "scans").exists(),
                   help="After changing targets or m/z windows"):
        run_steps(root, [("Extract traces", cli.cmd_extract, none), ("Analyze", cli.cmd_analyze, none)])
    if b[3].button("Re-analyze only", disabled=not status["Traces"],
                   help="After changing time windows or thresholds"):
        run_steps(root, [("Analyze", cli.cmd_analyze, none)])
    if st.session_state.pop("refresh", False):
        st.rerun()  # update the status list and the results tab
    if st.session_state.get("log"):
        st.code(st.session_state["log"])


# ------------------------------------------------------------------ results
def load(name):
    p = rep / name
    return pd.read_csv(p) if p.exists() else pd.DataFrame()


with tabs[5]:
    gw = load("global_windows.csv")
    res_all = load("trace_results.csv")
    if res_all.empty:
        st.info("No results yet. Run the analysis in tab 5.")
        st.stop()
    if gw.empty:
        gw = pd.DataFrame(columns=["ppm", "axis", "start", "length", "pct_all", "pct_with_stable", "n_traces"])
    c = st.columns([1, 1, 3])
    ppm_sel = c[0].selectbox("m/z window (± ppm)", sorted(res_all["ppm"].unique(), reverse=True),
                             format_func=lambda v: f"{v:g}")
    axis_sel = c[1].selectbox("Time axis", ["clock", "onset"],
                              format_func={"clock": "Acquisition time", "onset": "Time after signal onset"}.get)
    tag = f"ppm{ppm_sel:g}"
    res = res_all[(res_all["ppm"] == ppm_sel) & (res_all["axis"] == axis_sel)]
    g = gw[(gw["ppm"] == ppm_sel) & (gw["axis"] == axis_sel)]
    thr_all = json.loads((rep / "thresholds.json").read_text())
    thr = thr_all.get(tag, {})

    if thr.get("n_traces", 99) < 5:
        st.warning(f"Only {thr.get('n_traces')} trace(s) could be analyzed with this m/z window, so the thresholds "
                   "derived from the data are not reliable. See 'Detection and flags' below for the ions that "
                   "were not detected; add acquisitions, or set the thresholds manually.")
    st.subheader("Global window")
    m = st.columns(4)
    if len(g):
        r = g.iloc[0]
        m[0].metric("Recommended window", f"{r['start']:g} – {r['start'] + r['length']:g} s")
        m[1].metric("Traces stable in it", f"{r['pct_all']:.0f} %", help=f"of {int(r['n_traces'])} analyzed traces")
        m[2].metric("Among traces with a stable window", f"{r['pct_with_stable']:.0f} %")
    else:
        m[0].metric("Recommended window", "none")
    m[3].metric("Traces without stable window", f"{int(res['no_stable_window'].sum())} / {len(res)}")

    f = st.columns(2)
    for col, name in zip(f, (f"profiles_{tag}.png", f"stable_pct_{axis_sel}_{tag}.png")):
        if (rep / name).exists():
            col.image(str(rep / name))

    st.subheader("Thresholds derived from the data")
    st.caption("Inspect the figures before accepting the thresholds. To replace them, enter values and "
               "re-analyze; leave empty to use the derived values.")
    t1, t2, t3 = st.columns([1, 1, 2])
    ov = cfg["thresholds"].get("override") or {}
    t1.metric("RSD limit = noise model ×", f"{thr.get('rsd_factor', np.nan):.2f}")
    t2.metric("Slope limit, % per s", f"{thr.get('slope_limit', np.nan):.2f}")
    with t3.form("override"):
        o1, o2 = st.columns(2)
        new_f = o1.text_input("Override RSD factor", "" if ov.get("rsd_factor") is None else str(ov["rsd_factor"]))
        new_s = o2.text_input("Override slope limit", "" if ov.get("slope_limit") is None else str(ov["slope_limit"]))
        if st.form_submit_button("Apply and re-analyze"):
            try:
                cfg["thresholds"]["override"] = {"rsd_factor": opt_float(new_f), "slope_limit": opt_float(new_s)}
                save_cfg(root, cfg)
                if run_steps(root, [("Analyze", cli.cmd_analyze, argparse.Namespace())]):
                    st.rerun()
            except ValueError as e:
                st.error(f"Not applied: {e}")
    f = st.columns(2)
    for col, name in zip(f, (f"rsd_vs_intensity_{tag}.png", f"threshold_histograms_{tag}.png")):
        if (rep / name).exists():
            col.image(str(rep / name))

    st.subheader("Manual window")
    st.caption("Optional. Set a window on acquisition time (start and end, in s) for all ions, or for single "
               "ions in the table; a window set for one ion takes precedence. Manual windows are evaluated "
               "with the same criteria as the automatic ones and exported with them. Empty = no manual window.")
    manual = cfg.get("manual_windows") or {"global": None, "targets": {}}
    all_targets = sorted(load("trace_results.csv")["target_id"].unique())
    with st.form("manual"):
        g1, g2, _ = st.columns([1, 1, 2])
        gman = manual.get("global") or ["", ""]
        g_start = g1.text_input("All ions: start, s", str(gman[0]))
        g_end = g2.text_input("All ions: end, s", str(gman[1]))
        per = manual.get("targets") or {}
        man_df = pd.DataFrame({
            "target ion": all_targets,
            "start, s": [per.get(t, [None, None])[0] for t in all_targets],
            "end, s": [per.get(t, [None, None])[1] for t in all_targets],
        }).astype({"start, s": float, "end, s": float})
        with st.expander("Per ion (optional)", expanded=bool(per)):
            man_edit = st.data_editor(man_df, hide_index=True, width="stretch", disabled=["target ion"],
                                      key="manual_editor")
        if st.form_submit_button("Apply and re-analyze"):
            try:
                def pair(a, b, name):
                    if a is None or b is None or (isinstance(a, float) and np.isnan(a)) or (isinstance(b, float) and np.isnan(b)):
                        if not (a is None or (isinstance(a, float) and np.isnan(a))) or \
                                not (b is None or (isinstance(b, float) and np.isnan(b))):
                            raise ValueError(f"{name}: give both start and end")
                        return None
                    if not 0 <= a < b:
                        raise ValueError(f"{name}: start must be 0 or more and smaller than end")
                    return [float(a), float(b)]
                new_manual = {"global": pair(opt_float(g_start), opt_float(g_end), "All ions"), "targets": {}}
                for r in man_edit.to_dict("records"):
                    w = pair(r["start, s"], r["end, s"], r["target ion"])
                    if w:
                        new_manual["targets"][r["target ion"]] = w
                cfg["manual_windows"] = new_manual
                save_cfg(root, cfg)
                if run_steps(root, [("Analyze", cli.cmd_analyze, argparse.Namespace())]):
                    st.rerun()
            except ValueError as e:
                st.error(f"Not applied: {e}")
    choice = load("window_choice.csv")
    if len(choice):
        st.caption("Automatic and manual windows per ion, on acquisition time. pct_stable: acquisitions in "
                   "which the signal is stable in that window. Also saved as window_choice.csv.")
        ch = choice[choice["ppm"] == ppm_sel].drop(columns="ppm")
        st.dataframe(ch, hide_index=True, width="stretch",
                     column_config={"pct_stable": st.column_config.NumberColumn(format="%.0f"),
                                    "median_rsd_pct": st.column_config.NumberColumn(format="%.2f"),
                                    "median_abs_slope": st.column_config.NumberColumn(format="%.3f"),
                                    "mean_intensity": st.column_config.NumberColumn(format="%.3e"),
                                    "cv_pct_between_acquisitions": st.column_config.NumberColumn(format="%.2f")})
        st.download_button("Download windows (.csv)", choice.to_csv(index=False), "window_choice.csv", "text/csv")

    st.subheader("Per target ion")
    grp = load("group_windows.csv")
    grp = grp[(grp["ppm"] == ppm_sel) & (grp["axis"] == axis_sel)].drop(columns=["ppm", "axis"])
    st.dataframe(grp, hide_index=True, width="stretch")
    targets_av = sorted(res["target_id"].unique())
    if targets_av and (work / "traces.parquet").exists():
        tid = st.selectbox("Target ion", targets_av)
        tr = pd.read_parquet(work / "traces.parquet", filters=[("ppm", "==", ppm_sel), ("target_id", "==", tid)])
        tr = tr[tr["file"].isin(res.loc[res["target_id"] == tid, "file"])]
        sel = res[(res["target_id"] == tid) & ~res["no_stable_window"]][["file", "start", "length", "onset_s"]]
        origin = sel["onset_s"] if axis_sel == "onset" else 0.0
        sel = sel.assign(t0=sel["start"] + origin, t1=sel["start"] + sel["length"] + origin)
        inwin = tr.merge(sel, on="file")
        inwin = inwin[(inwin["time_s"] >= inwin["t0"]) & (inwin["time_s"] < inwin["t1"])]
        x = alt.X("time_s:Q", title="Acquisition time (s)")
        y = alt.Y("intensity:Q", title="Absolute intensity", scale=alt.Scale(zero=True))
        tip = [alt.Tooltip("file:N"), alt.Tooltip("time_s:Q", format=".1f"), alt.Tooltip("intensity:Q", format=".3s")]
        base = alt.Chart(tr.dropna(subset=["intensity"])).mark_line(color=GREY, strokeWidth=1, opacity=0.8).encode(
            x=x, y=y, detail="file:N", tooltip=tip)
        hi = alt.Chart(inwin.dropna(subset=["intensity"])).mark_line(color=BLUE, strokeWidth=2).encode(
            x=x, y=y, detail="file:N", tooltip=tip)
        chart = base + hi
        mw = manual_window_for(tid, manual)
        if mw:
            rules = alt.Chart(pd.DataFrame({"x": [mw[0], mw[1]]})).mark_rule(
                color="#eb6834", strokeWidth=2, strokeDash=[5, 4]).encode(x="x:Q")
            chart = chart + rules
        st.altair_chart(chart.properties(height=340).interactive(bind_y=False), width="stretch")
        st.caption("One line per acquisition. Blue: the stable window selected for that acquisition. "
                   "Acquisitions with no blue segment have no stable window."
                   + (f" Dashed orange lines: manual window ({mw[0]:g}–{mw[1]:g} s, set for {mw[2]})." if mw else ""))
        show = ["file", "start", "length", "mean", "rsd", "rel_slope", "no_stable_window", "plateau_start_s",
                "plateau_end_s", "trend_rel_slope", "ppm_error_median", "iso_status", "saturated", "it_at_max",
                "unresolved", "multiple_centroids"]
        st.dataframe(res.loc[res["target_id"] == tid, [c for c in show if c in res]], hide_index=True,
                     width="stretch")

    st.subheader("Detection and flags")
    summ = load("trace_summary.csv")
    summ = summ[summ["ppm"] == ppm_sel]
    flags = [c for c in ("out_of_range", "not_detected", "isotope_mismatch", "saturated", "it_at_max",
                         "unresolved", "multiple_centroids", "analyzed") if c in summ]
    st.dataframe(summ.groupby("target_id")[flags].agg(lambda s: int(s.fillna(False).astype(bool).sum()))
                 .assign(acquisitions=summ.groupby("target_id").size()).reset_index(),
                 hide_index=True, width="stretch")

    cmp_ = load("window_comparison.csv")
    if len(cmp_):
        st.subheader("Comparison of m/z windows")
        st.caption("lost_fraction: scans detected with the wide window but not the narrow one. "
                   "interference_removed_fraction: scans with several peaks in the wide window and one in the narrow.")
        a, bcol = st.columns([1, 2])
        if (rep / "mass_error.png").exists():
            a.image(str(rep / "mass_error.png"))
        bcol.dataframe(cmp_, hide_index=True, width="stretch")

    val = load("validation.csv")
    if len(val):
        st.subheader("Validation (windows derived on half of the acquisitions, tested on the rest)")
        st.dataframe(val[(val["ppm"] == ppm_sel) & (val["axis"] == axis_sel)].drop(columns=["ppm", "axis"]),
                     hide_index=True, width="stretch")

    with st.expander("Dependence on intensity, m/z, reaction and run order"):
        mm = rep / f"mixed_model_{tag}.txt"
        st.code(mm.read_text() if mm.exists() else "not available")

    if st.button("Prepare report archive"):
        st.session_state["zip"] = Path(shutil.make_archive(str(work / "report"), "zip", rep)).read_bytes()
    if "zip" in st.session_state:
        st.download_button("Download report (.zip)", st.session_state["zip"], "mzstab_report.zip")
