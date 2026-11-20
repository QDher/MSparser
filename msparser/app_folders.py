"""Interface of the intensity extraction, organised by data folder.

The user gives a folder (with or without its subfolders). Every folder that
holds acquisitions gets an `MSparser` subfolder with its converted files,
compound list, ion list, settings and results (see folders.py).
"""
from __future__ import annotations

import shutil
from pathlib import Path

import altair as alt
import numpy as np
import pandas as pd
import streamlit as st

from msparser import cli
from msparser import folders as fo
from msparser import ui_state
from msparser.intensity import FLAGS
from msparser.ui_shared import components_editor, converter_panel, read_components, targets_editor

GREY, BLUE = "#898781", "#2a78d6"
PEAK_COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
FLAG_HELP = {
    "low_stability": "RSD in the window above the limit",
    "window_variation": "the time windows of the acquisition differ by more than the limit",
    "gaps": "signal missing in more than 10 % of the scans of the window",
    "multiple_peaks": "several peaks in the m/z window",
    "wide_mz_window": "found only with an m/z window wider than the first one",
    "isotope_mismatch": "isotope peak ratio differs from the theoretical one",
    "saturated": "intensity at the saturation level",
    "short_acquisition": "too few scans per window; all scans of the acquisition used",
}
VIEWS = {"Intensity": ("", "%.3e"), "m/z detected": ("_mz_detected", "%.5f"), "Mass error (ppm)": ("_ppm_error", "%.2f")}
NUM = {"intensity": st.column_config.NumberColumn(format="%.3e"),
       "mz_theoretical": st.column_config.NumberColumn("m/z theoretical", format="%.5f"),
       "rsd_pct": st.column_config.NumberColumn("RSD %", format="%.1f"),
       "window_variation_pct": st.column_config.NumberColumn("variation %", format="%.0f"),
       "ppm_error": st.column_config.NumberColumn(format="%.2f"),
       "mz_observed": st.column_config.NumberColumn("m/z detected", format="%.5f")}


def _load(path: Path) -> pd.DataFrame:
    return pd.read_csv(path) if Path(path).exists() else pd.DataFrame()


def _parse_numbers(text: str) -> list[float]:
    vals = [float(x) for x in text.replace(";", ",").split(",") if x.strip()]
    if not vals or min(vals) <= 0:
        raise ValueError("enter positive numbers separated by commas")
    return vals


# ================================================================== entry
def render(arg_folder: str | None, state: dict) -> None:
    # ------------------------------------------------------------ sidebar
    if "data_folder_input" not in st.session_state:
        st.session_state["data_folder_input"] = arg_folder or state.get("data_folder") or ""
    typed = st.sidebar.text_input(
        "Data folder", key="data_folder_input", placeholder=r"C:\data\experiment1",
        help="The folder that holds your acquisitions (raw, mzML or CSV files). Everything the program "
             "writes goes into an 'MSparser' subfolder of each data folder.")
    include = st.sidebar.checkbox("Include subfolders", value=bool(state.get("include_subfolders", True)))
    recent = [f for f in state.get("recent_folders", []) if f != typed.strip().strip('"')]
    if recent:
        def _pick():
            if st.session_state.get("recent_pick"):
                st.session_state["data_folder_input"] = st.session_state["recent_pick"]
                st.session_state["recent_pick"] = None
        st.sidebar.selectbox("Recent folders", recent, index=None, key="recent_pick", on_change=_pick,
                             placeholder="Choose a recent folder")
    folder = typed.strip().strip('"')
    if not folder:
        st.info("Enter in the left panel the folder that holds your acquisitions. Tick 'Include subfolders' to "
                "process its subfolders as well. Tip: in File Explorer, hold Shift, right-click the folder and "
                "choose 'Copy as path'.")
        return
    root = Path(folder).expanduser()
    if not root.is_dir():
        st.sidebar.warning("Folder not found.")
        st.warning(f"The folder `{root}` does not exist.")
        return
    root = root.resolve()
    if state.get("data_folder") != str(root) or state.get("include_subfolders") != include:
        state.update(data_folder=str(root), include_subfolders=include)
        ui_state.remember_folder(state, str(root))
        ui_state.save(state)

    settings = fo.load_settings(root)
    ic = settings["intensity"]

    def save():
        fo.save_settings(root, settings)

    found = fo.find_units(root, include)
    rel = {u: u.relative_to(root).as_posix() for u in found}
    excluded = set(settings.get("excluded_folders") or [])
    units = [u for u in found if rel[u] not in excluded]
    labels = {u: fo.label(u, root) for u in found}
    stat = {u: fo.unit_status(u) for u in found}
    st.sidebar.write(f"{'✅' if units else '⬜'} {len(units)} folder(s) used")
    st.sidebar.write(f"{'✅' if units and all(stat[u]['compounds'] for u in units) else '⬜'} "
                     f"{sum(bool(stat[u]['compounds']) for u in units)} with a compound list")
    st.sidebar.write(f"{'✅' if units and all(stat[u]['results'] for u in units) else '⬜'} "
                     f"{sum(stat[u]['results'] for u in units)} with results")

    tabs = st.tabs(["1. Folders", "2. Compounds", "3. Ions", "4. Convert and extract", "5. Results"])

    # ------------------------------------------------------------ folders
    with tabs[0]:
        st.caption("Each folder that holds acquisitions is processed on its own. The program writes into an "
                   f"`{fo.SUB}` subfolder of each one: the converted files, the compound list, the ions, the "
                   "settings used and the results. Your raw files are not modified.")
        if not found:
            st.warning("No raw, mzML or CSV file found in this folder"
                       + (" or its subfolders." if include else ". Tick 'Include subfolders' to look deeper."))
        else:
            table = pd.DataFrame({
                "use": [rel[u] not in excluded for u in found],
                "folder": [labels[u] for u in found],
                "raw files": [stat[u]["raw_files"] for u in found],
                "converted or mzML/CSV files": [stat[u]["usable_files"] for u in found],
                "compounds": [stat[u]["compounds"] for u in found],
                "results": ["yes" if stat[u]["results"] else "" for u in found],
            })
            edited = st.data_editor(table, hide_index=True, width="stretch", key=f"units_{root}",
                                    disabled=[c for c in table.columns if c != "use"])
            now_excluded = sorted(rel[u] for u, use in zip(found, edited["use"]) if not bool(use))
            if now_excluded != sorted(excluded):
                settings["excluded_folders"] = now_excluded
                save()
                st.rerun()
            st.caption(f"Folders are written from `{root.parent}`. Untick a folder to leave it out. "
                       f"{len(units)} of {len(found)} folder(s) used.")
        converter_panel(settings, save, "fold")

    if not units:
        for t in tabs[1:]:
            with t:
                st.info("No folder to work on: see tab 1.")
        return

    def pick(key: str, label: str = "Folder") -> Path:
        names = [labels[u] for u in units]
        chosen = st.selectbox(label, names, key=key) if len(units) > 1 else names[0]
        if len(units) == 1:
            st.caption(f"Folder: `{chosen}`")
        return units[names.index(chosen)]

    # ---------------------------------------------------------- compounds
    with tabs[1]:
        st.caption("Each folder has its own compound list. Enter the list of one folder, then copy it to the "
                   "others if they contain the same compounds.")
        u = pick("comp_folder")
        comp_path = u / fo.SUB / "components.yaml"
        components_editor(comp_path, settings["targets"]["pairs"], f"comp_{rel[u]}", "tab 3")
        others = [x for x in units if x != u]
        if others:
            st.divider()
            c1, c2 = st.columns([1, 2])
            if c1.button(f"Copy this list to the {len(others)} other folder(s)", disabled=not comp_path.exists(),
                         help="Replaces the compound list of every other folder used."):
                n = fo.copy_components(u, others)
                st.session_state["copied"] = f"Compound list of `{labels[u]}` copied to {n} folder(s)."
                st.rerun()
            if st.session_state.get("copied"):
                c2.success(st.session_state.pop("copied"))
            overview = pd.DataFrame({
                "folder": [labels[x] for x in units],
                "compounds": [", ".join(d["name"] for d in read_components(x / fo.SUB / "components.yaml")) or "(none)"
                              for x in units]})
            st.dataframe(overview, hide_index=True, width="stretch")

    # --------------------------------------------------------------- ions
    with tabs[2]:
        st.caption("The ions followed in each folder, computed from its compound list. They are generated "
                   "automatically at extraction; generate them here to review them and untick some.")
        u = pick("ion_folder")
        cfg_u = fo.unit_cfg(u, settings)
        targets_editor(cfg_u, f"ions_{rel[u]}")
        comp, work = Path(cfg_u["paths"]["components"]), Path(cfg_u["paths"]["work_dir"])
        if comp.exists() and (work / "targets.csv").exists():
            # keep the mark used at extraction in step, so that a reviewed selection is not regenerated
            (work / "components.sha1").write_text(fo._components_hash(comp))

    # ------------------------------------------------------------ extract
    with tabs[3]:
        st.caption("For each ion in each file: the m/z windows are tried in increasing order and the first one "
                   f"containing a signal is used. The acquisition is covered by {int(ic['n_windows'])} time "
                   f"windows, each lasting {ic['window_fraction']:.0%} of the acquisition and spread evenly from "
                   "start to end (for 30 s: six 10 s windows starting every 4 s). The intensity is the mean "
                   "over the most stable of them. Nothing is excluded for instability; such values are flagged.")
        with st.form("fold_settings"):
            c1, c2, c3 = st.columns(3)
            ppm = c1.text_input("m/z windows tried, ± ppm (in this order)", ", ".join(f"{p:g}" for p in ic["ppm_steps"]))
            options = sorted({1, 3, 5, int(ic["min_scans"])})
            min_scans = c1.selectbox("Consecutive scans needed to report a signal", options,
                                     index=options.index(int(ic["min_scans"])))
            variation = c2.number_input("Flag when windows differ by more than, %",
                                        value=float(ic["variation_flag_pct"]), min_value=1.0, step=5.0,
                                        help="Difference between the highest and lowest window mean of the "
                                             "acquisition, relative to the reported intensity.")
            rsd_flag = c2.number_input("Flag low stability above RSD, %", value=float(ic["rsd_flag_pct"]),
                                       min_value=1.0, step=5.0)
            rel_min = c3.number_input("Report other peaks from, % of the most intense",
                                      value=float(ic["min_relative_intensity"]) * 100, min_value=1.0,
                                      max_value=100.0, step=5.0)
            if st.form_submit_button("Save settings"):
                try:
                    ic.update(ppm_steps=sorted(_parse_numbers(ppm)), min_scans=int(min_scans),
                              variation_flag_pct=float(variation), rsd_flag_pct=float(rsd_flag),
                              min_relative_intensity=float(rel_min) / 100)
                    save()
                    st.success("Settings saved. They apply to all folders.")
                except ValueError as e:
                    st.error(f"Not saved: {e}")
        missing = [labels[x] for x in units if not stat[x]["compounds"]]
        if missing:
            st.warning("No compound list in: " + ", ".join(missing) + ". These folders will be skipped (tab 2).")
        b1, b2, _ = st.columns([1, 1, 2])
        go = b1.button("Convert and extract", type="primary",
                       help="Converts the raw files not yet converted, then extracts the intensities. Files "
                            "already processed with the same settings and ions are reused.")
        go_no_conv = b2.button("Extract only", help="Uses the files already converted; no conversion.")
        if go or go_no_conv:
            bar = st.progress(0.0, text="Starting...")

            def progress(i, n, f, j, m, name):
                frac = ((i - 1) + (j / m if m else 0)) / n
                bar.progress(min(frac, 1.0), text=f"Folder {i} / {n}: {fo.label(f, root)}"
                             + (f" — file {j} / {m}" if m else " — converting"))

            res = fo.run_all(root, units, settings, convert=go, progress=progress)
            bar.empty()
            st.session_state["fold_run"] = res.to_dict("records")
            st.rerun()
        if st.session_state.get("fold_run"):
            res = pd.DataFrame(st.session_state["fold_run"])
            ok = int((res["status"] == "done").sum())
            (st.success if ok == len(res) else st.warning)(f"{ok} of {len(res)} folder(s) done.")
            st.dataframe(res[["data_folder", "status", "message", "files", "values", "flagged", "not_detected"]],
                         hide_index=True, width="stretch")
            with st.expander("Log"):
                st.code("\n".join(f"== {r['data_folder']} ==\n{r['log']}" for r in res.to_dict("records")))

    # ------------------------------------------------------------ results
    with tabs[4]:
        done = [x for x in units if stat[x]["results"]]
        if not done:
            st.info("No results yet. Convert and extract in tab 4.")
            return
        scopes = (["All folders"] if len(done) > 1 else []) + [labels[x] for x in done]
        scope = st.selectbox("Results of", scopes, key="res_scope") if len(scopes) > 1 else scopes[0]
        if scope == "All folders":
            _results_all(root, done)
        else:
            u = done[[labels[x] for x in done].index(scope)]
            _results_unit(u, labels[u], settings, save, root, done, rel[u])


# ===================================================== results, all folders
def _results_all(root: Path, done: list[Path]) -> None:
    out = root / fo.SUB / "results_all"
    if not (out / "intensity_values_all.csv").exists():
        fo.combine(root, done)
    values, matrix = _load(out / "intensity_values_all.csv"), _load(out / "intensity_matrix_all.csv")
    if values.empty:
        st.info("No combined results yet.")
        return
    values["flags"] = values["flags"].fillna("")
    det = values[values["status"] == "detected"]
    m = st.columns(5)
    m[0].metric("Folders", values["data_folder"].nunique())
    m[1].metric("Files", len(values[["data_folder", "file"]].drop_duplicates()))
    m[2].metric("Ions detected", det["target_id"].nunique(), help="in at least one file")
    m[3].metric("Intensity values", len(det))
    m[4].metric("Flagged values", int((det["flags"] != "").sum()))
    st.caption(f"Combined tables, saved in `{out}`. Choose one folder above for the choice of peaks and the "
               "detail of one ion in one file.")
    fl = st.columns([2, 2, 1])
    folders = sorted(values["data_folder"].unique())
    pick_f = fl[0].multiselect("Folders", folders, default=folders, key="all_folders")
    ions = sorted(det["target_id"].unique())
    pick_i = fl[1].multiselect("Ions", ions, default=ions, key="all_ions")
    layout = fl[2].radio("Layout", ["files × ions", "one row per value"], key="all_layout")
    if layout == "files × ions":
        view = st.radio("Values shown", list(VIEWS), horizontal=True, key="all_view")
        suffix, fmt = VIEWS[view]
        matrix = _load(out / f"intensity_matrix{suffix}_all.csv") if suffix else matrix
        mx = matrix[matrix["data_folder"].isin(pick_f)][["data_folder", "file"] + [c for c in pick_i if c in matrix]]
        st.dataframe(mx, hide_index=True, width="stretch",
                     column_config={c: st.column_config.NumberColumn(format=fmt) for c in pick_i})
        st.caption("Empty cell: ion not detected in that file, or not in the ion list of that folder.")
        st.download_button("Download this table (.csv)", mx.to_csv(index=False),
                           f"intensity_matrix{suffix}_all.csv", "text/csv")
    else:
        only = st.checkbox("Only flagged values", key="all_only")
        lv = values[values["data_folder"].isin(pick_f)
                    & (values["target_id"].isin(pick_i) | (values["status"] != "detected"))]
        if only:
            lv = lv[lv["flags"] != ""]
        show = ["data_folder", "file", "target_id", "status", "intensity", "mz_theoretical", "mz_observed",
                "ppm_error", "flags", "rsd_pct", "window_variation_pct", "mz_window_ppm", "window_start_s",
                "window_end_s", "selected_by"]
        st.dataframe(lv[[c for c in show if c in lv]], hide_index=True, width="stretch", column_config=NUM)
        st.download_button("Download this table (.csv)", lv.to_csv(index=False), "intensity_values_all.csv", "text/csv")
    with st.expander("Meaning of the flags"):
        st.markdown("\n".join(f"- **{k}**: {v}" for k, v in FLAG_HELP.items() if k in FLAGS))


# ======================================================= results, one folder
def _results_unit(u: Path, name: str, settings: dict, save, root: Path, done: list[Path], key: str) -> None:
    ic = settings["intensity"]
    cfg = fo.unit_cfg(u, settings)
    rep, work = Path(cfg["paths"]["report_dir"]), Path(cfg["paths"]["work_dir"])
    values = _load(rep / "intensity_values.csv")
    files, peaks, matrix = (_load(rep / f) for f in ("intensity_files.csv", "intensity_peaks.csv", "intensity_matrix.csv"))
    values["flags"] = values["flags"].fillna("")
    det = values[values["status"] == "detected"]
    failed = files[files["processing"] == "failed"] if len(files) else files
    m = st.columns(5)
    m[0].metric("Files", len(files), help=f"{len(failed)} could not be read")
    m[1].metric("Ions detected", det["target_id"].nunique(), help="in at least one file")
    m[2].metric("Intensity values", len(det))
    m[3].metric("Flagged values", int((det["flags"] != "").sum()))
    m[4].metric("Not detected", int((values["status"] == "not detected").sum()),
                help="ion and file combinations with no signal in any m/z window")
    st.caption(f"Tables of this folder are saved in `{rep}`.")
    if len(failed):
        with st.expander(f"{len(failed)} file(s) could not be read"):
            st.dataframe(failed[["file", "error"]], hide_index=True, width="stretch")

    st.subheader("Intensity table")
    fl = st.columns([3, 1])
    ions = sorted(det["target_id"].unique())
    pick_i = fl[0].multiselect("Ions", ions, default=ions, key=f"ions_{key}")
    layout = fl[1].radio("Layout", ["files × ions", "one row per value"], key=f"layout_{key}")
    if layout == "files × ions":
        view = st.radio("Values shown", list(VIEWS), horizontal=True, key=f"view_{key}")
        suffix, fmt = VIEWS[view]
        matrix = _load(rep / f"intensity_matrix{suffix}.csv") if suffix else matrix
        mx = matrix[["file"] + [c for c in pick_i if c in matrix]] if len(matrix) else matrix
        st.dataframe(mx, hide_index=True, width="stretch",
                     column_config={c: st.column_config.NumberColumn(format=fmt) for c in pick_i})
        st.caption("Empty cell: ion not detected in that file. Flags are in the 'one row per value' layout.")
        st.download_button("Download this table (.csv)", mx.to_csv(index=False),
                           f"intensity_matrix{suffix}.csv", "text/csv", key=f"dl_m_{key}")
    else:
        only = st.checkbox("Only flagged values", key=f"only_{key}")
        lv = values[values["target_id"].isin(pick_i) | (values["status"] != "detected")]
        if only:
            lv = lv[lv["flags"] != ""]
        show = ["file", "target_id", "status", "intensity", "mz_theoretical", "mz_observed", "ppm_error", "flags",
                "rsd_pct", "window_variation_pct", "mz_window_ppm", "window_start_s", "window_end_s",
                "scans_in_window", "selected_by"]
        st.dataframe(lv[[c for c in show if c in lv]], hide_index=True, width="stretch", column_config=NUM)
        st.download_button("Download this table (.csv)", lv.to_csv(index=False), "intensity_values.csv", "text/csv",
                           key=f"dl_v_{key}")
    with st.expander("Meaning of the flags"):
        st.markdown("\n".join(f"- **{k}**: {v}" for k, v in FLAG_HELP.items() if k in FLAGS))
    ion_tab = _load(rep / "intensity_ions.csv")
    if len(ion_tab):
        with st.expander("Per ion: theoretical and detected m/z, mass error"):
            st.caption("m/z detected is the median m/z of the peak over its scans; the columns below summarise it "
                       "over the files of this folder. Saved as intensity_ions.csv.")
            st.dataframe(ion_tab[ion_tab["files_detected"] > 0], hide_index=True, width="stretch", column_config={
                "mz_theoretical": st.column_config.NumberColumn("m/z theoretical", format="%.5f"),
                "mz_observed_median": st.column_config.NumberColumn("m/z detected (median)", format="%.5f"),
                "ppm_error_median": st.column_config.NumberColumn("ppm error (median)", format="%.2f"),
                "ppm_error_min": st.column_config.NumberColumn("ppm min", format="%.2f"),
                "ppm_error_max": st.column_config.NumberColumn("ppm max", format="%.2f"),
                "intensity_median": st.column_config.NumberColumn("intensity (median)", format="%.3e")})

    # ----------------------------------------------------- choice of the peak
    st.subheader("Choice of the peak when several are found")
    st.caption("When the m/z window holds several peaks, the one closest to the theoretical m/z, the most "
               "intense one and those of at least the set percentage of the most intense are listed. "
               "Tick the peak to export; unticked groups follow the rule below, which applies to all folders.")
    rules = {"closest": "Closest to the theoretical m/z", "most_intense": "Most intense"}
    rule = st.radio("Rule for all other ions", list(rules), index=list(rules).index(ic["default_peak"]),
                    format_func=rules.get, horizontal=True, key=f"rule_{key}")
    multi = peaks[(peaks["status"] == "detected") & (peaks["n_peaks_in_window"] > 1)].copy() if len(peaks) else peaks
    edited = None
    if len(multi):
        multi = multi.sort_values(["file", "target_id", "mz_observed"])
        cols = ["selected", "file", "target_id", "mz_observed", "ppm_error", "intensity",
                "relative_to_most_intense_pct", "rsd_pct", "is_closest", "is_most_intense"]
        edited = st.data_editor(
            multi[cols].rename(columns={"selected": "use"}), hide_index=True, width="stretch",
            disabled=[c for c in cols if c != "selected"], key=f"peak_editor_{key}",
            column_config={"mz_observed": st.column_config.NumberColumn("m/z", format="%.5f"),
                           "ppm_error": st.column_config.NumberColumn("ppm", format="%.2f"),
                           "intensity": st.column_config.NumberColumn(format="%.3e"),
                           "relative_to_most_intense_pct": st.column_config.NumberColumn("% of most intense", format="%.0f"),
                           "rsd_pct": st.column_config.NumberColumn("RSD %", format="%.1f"),
                           "is_closest": "closest", "is_most_intense": "most intense"})
    else:
        st.info("No ion has more than one peak in its m/z window.")
    if st.button("Apply the choice of peaks", key=f"apply_{key}"):
        try:
            rows = []
            if edited is not None:
                flag = "is_most_intense" if rule == "most_intense" else "is_closest"
                for (file, tid), g in edited.groupby(["file", "target_id"], sort=False):
                    used = g[g["use"].astype(bool)]
                    if len(used) > 1:
                        raise ValueError(f"{file}, {tid}: tick one peak only")
                    if len(used) == 1 and not bool(used[flag].iloc[0]):  # differs from the rule
                        rows.append(dict(file=file, target_id=tid, mz_observed=float(used["mz_observed"].iloc[0])))
            sel_path = work / "peak_selection.csv"
            if rows:
                pd.DataFrame(rows).to_csv(sel_path, index=False)
            elif sel_path.exists():
                sel_path.unlink()
            ic["default_peak"] = rule
            save()
            for x in done:  # the rule is common to all folders
                cli.write_intensity_outputs(fo.unit_cfg(x, settings, write=True))
            fo.combine(root, done)
            st.session_state.pop(f"peak_editor_{key}", None)
            st.rerun()
        except ValueError as e:
            st.error(f"Not applied: {e}")

    # ---------------------------------------------------------------- detail
    st.subheader("Detail of one ion in one file")
    det_files = sorted(det["file"].unique())
    if not det_files:
        return
    d1, d2 = st.columns(2)
    f_sel = d1.selectbox("File", det_files, key=f"file_{key}")
    t_sel = d2.selectbox("Ion", sorted(det.loc[det["file"] == f_sel, "target_id"]), key=f"ion_{key}")
    row = det[(det["file"] == f_sel) & (det["target_id"] == t_sel)].iloc[0]
    k = st.columns(5)
    k[0].metric("Intensity", f"{row['intensity']:.3e}")
    k[1].metric("RSD in window", f"{row['rsd_pct']:.1f} %")
    k[2].metric("Variation between windows", f"{row['window_variation_pct']:.0f} %",
                help="(highest − lowest window mean) / reported intensity")
    k[3].metric("m/z window used", f"± {row['mz_window_ppm']:g} ppm")
    k[4].metric("Mass error", f"{row['ppm_error']:.1f} ppm")
    if row["flags"]:
        st.warning("Flags: " + "; ".join(f"{f} ({FLAG_HELP.get(f, '')})" for f in row["flags"].split("; ")))
    g1, g2 = st.columns(2)
    tr_path = work / "intensity_traces.parquet"
    if tr_path.exists():
        tr = pd.read_parquet(tr_path, filters=[("file", "==", f_sel), ("target_id", "==", t_sel)])
        if len(tr):
            order = sorted(tr["mz_observed"].unique())
            tr["peak"] = tr["mz_observed"].map(lambda v: f"m/z {v:.5f}")
            names = [f"m/z {v:.5f}" for v in order]
            band = alt.Chart(pd.DataFrame({"x": [row["window_start_s"]], "x2": [row["window_end_s"]]})).mark_rect(
                color=BLUE, opacity=0.10).encode(x="x:Q", x2="x2:Q")
            enc = dict(x=alt.X("time_s:Q", title="Acquisition time (s)"),
                       y=alt.Y("intensity:Q", title="Absolute intensity", scale=alt.Scale(zero=True)),
                       tooltip=[alt.Tooltip("peak:N"), alt.Tooltip("time_s:Q", format=".1f"),
                                alt.Tooltip("intensity:Q", format=".3s")])
            if len(order) > 1:
                enc["color"] = alt.Color("peak:N", title="Peak", legend=alt.Legend(orient="bottom"),
                                         scale=alt.Scale(domain=names, range=PEAK_COLORS[:len(names)]))
                line = alt.Chart(tr).mark_line(strokeWidth=1.6).encode(**enc)
            else:
                line = alt.Chart(tr).mark_line(strokeWidth=1.6, color=BLUE).encode(**enc)
            with g1:
                st.markdown("**Signal over the acquisition**")
                st.altair_chart((band + line).properties(height=300), width="stretch")
                st.caption(f"Shaded: the window used ({row['window_start_s']:g}–{row['window_end_s']:g} s), "
                           f"the most stable of the windows. Exported peak: m/z {row['mz_observed']:.5f}.")
    wm = _load(rep / "intensity_window_means.csv")
    if len(wm):
        w = wm[(wm["file"] == f_sel) & (wm["target_id"] == t_sel)
               & np.isclose(wm["mz_observed"], row["mz_observed"], rtol=0, atol=1e-6)].copy()
        if len(w):
            w["mid"] = (w["start"] + w["end"]) / 2
            w["rsd_pct"] = w["rsd"] * 100
            w["window"] = np.where(w["selected"], "Window used", "Other windows")
            pts = alt.Chart(w).mark_circle(size=70).encode(
                x=alt.X("mid:Q", title="Centre of the window (s)"),
                y=alt.Y("mean:Q", title="Window mean intensity", scale=alt.Scale(zero=True)),
                color=alt.Color("window:N", title=None, legend=alt.Legend(orient="bottom"),
                                scale=alt.Scale(domain=["Window used", "Other windows"], range=[BLUE, GREY])),
                tooltip=[alt.Tooltip("start:Q", title="start (s)"), alt.Tooltip("end:Q", title="end (s)"),
                         alt.Tooltip("mean:Q", format=".3s"), alt.Tooltip("rsd_pct:Q", title="RSD %", format=".1f"),
                         alt.Tooltip("n:Q", title="scans")])
            with g2:
                st.markdown("**Variation between windows**")
                st.altair_chart(pts.properties(height=300), width="stretch")
                st.caption(f"One dot per window ({len(w)} windows of {row['window_end_s'] - row['window_start_s']:.1f} s). "
                           f"Lowest window mean {w['mean'].min():.3e}, highest {w['mean'].max():.3e}: a spread of "
                           f"{row['window_variation_pct']:.0f} % of the reported intensity "
                           f"(flagged from {ic['variation_flag_pct']:g} %).")
    st.divider()
    if st.button("Prepare archive of the tables of this folder", key=f"zip_{key}"):
        st.session_state[f"zipdata_{key}"] = Path(shutil.make_archive(str(work / "intensity_report"), "zip", rep)).read_bytes()
    if f"zipdata_{key}" in st.session_state:
        st.download_button("Download (.zip)", st.session_state[f"zipdata_{key}"], f"MSparser_{u.name}.zip",
                           key=f"zipdl_{key}")
