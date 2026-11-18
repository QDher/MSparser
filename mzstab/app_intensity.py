"""Interface of the intensity-extraction mode (called from app.py)."""
from __future__ import annotations

import contextlib
import io
import shutil
from pathlib import Path

import altair as alt
import numpy as np
import pandas as pd
import streamlit as st

from mzstab import cli
from mzstab.convert import find_raw
from mzstab.intensity import FLAGS, discover_files

GREY, BLUE, ORANGE = "#898781", "#2a78d6", "#eb6834"
PEAK_COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
FLAG_HELP = {
    "low_stability": "RSD in the window above the limit",
    "window_variation": "the 10 s windows of the acquisition differ by more than the limit",
    "gaps": "signal missing in more than 10 % of the scans of the window",
    "multiple_peaks": "several peaks in the m/z window",
    "wide_mz_window": "found only with an m/z window wider than the first one",
    "isotope_mismatch": "isotope peak ratio differs from the theoretical one",
    "saturated": "intensity at the saturation level",
    "short_acquisition": "acquisition shorter than the window; all scans used",
}


def _load(path: Path) -> pd.DataFrame:
    return pd.read_csv(path) if path.exists() else pd.DataFrame()


def render(ctx: dict):
    root, cfg, paths = ctx["root"], ctx["cfg"], ctx["paths"]
    work, rep = ctx["work"], ctx["rep"]
    save_cfg, resolved, parse_numbers = ctx["save_cfg"], ctx["resolved"], ctx["parse_numbers"]
    ic = cfg["intensity"]

    status = {
        "Components": Path(paths["components"]).exists(),
        "Targets": (work / "targets.csv").exists(),
        "Intensities": (rep / "intensity_values.csv").exists(),
    }
    for name, done in status.items():
        st.sidebar.write(("✅ " if done else "⬜ ") + name)

    tabs = st.tabs(["1. Components", "2. Files", "3. Targets", "4. Extract", "5. Results"])

    with tabs[0]:
        ctx["components_editor"]("tab 3")

    # ------------------------------------------------------------------ files
    with tabs[1]:
        st.caption("All files below these folders are used, including those in subfolders, so the data can be "
                   "organised in any folder tree. The same components apply to every file; each file only "
                   "reports the ions matching its own polarity.")
        st.markdown(f"Give a subfolder of the project folder `{root}`, or a full path.")
        f1, f2 = st.columns(2)
        in_raw = f1.text_input("Folder containing the raw files (with subfolders)", cfg["paths"]["raw_dir"],
                               key="int_raw")
        in_mzml = f2.text_input("Folder for the converted files (mzML), or containing mzML/CSV files",
                                cfg["paths"]["mzml_dir"], key="int_mzml",
                                help="The subfolders of the raw folder are reproduced here at conversion.")
        chosen = {"raw_dir": in_raw.strip().strip('"') or ".", "mzml_dir": in_mzml.strip().strip('"') or "."}
        if any(chosen[k] != cfg["paths"][k] for k in chosen):
            cfg["paths"].update(chosen)
            save_cfg(root, cfg)
            st.rerun()
        raw_p, mzml_p = Path(paths["raw_dir"]), Path(paths["mzml_dir"])
        n_raw = len(find_raw(raw_p, recursive=True)) if raw_p.is_dir() else 0
        (f1.caption if raw_p.is_dir() else f1.warning)(
            f"{raw_p} — {n_raw} raw file(s) found" if raw_p.is_dir() else f"{raw_p} — folder not found")
        found = discover_files(mzml_p)
        f2.caption(f"{mzml_p} — {len(found)} mzML/CSV file(s) found")
        if st.button("Convert raw files", disabled=n_raw == 0,
                     help="Needs ThermoRawFileParser or msconvert (program path under Stability analysis > "
                          "Settings > Advanced settings > Converters). Already converted files are skipped."):
            buf = io.StringIO()
            try:
                with st.spinner("Converting..."), contextlib.redirect_stdout(buf):
                    cli.cmd_convert(resolved(root), None)
                st.success("Conversion done.")
            except (Exception, SystemExit) as e:
                st.error(f"Conversion failed: {e}")
            st.session_state["int_log"] = buf.getvalue()
            found = discover_files(mzml_p)
        if st.session_state.get("int_log"):
            st.code(st.session_state["int_log"])
        if len(found):
            per_folder = found.groupby("folder").size().rename("files").reset_index()
            per_folder["folder"] = per_folder["folder"].replace("", "(top folder)")
            c1, c2 = st.columns([1, 2])
            c1.dataframe(per_folder, hide_index=True, width="stretch")
            c2.dataframe(found[["folder", "name"]], hide_index=True, width="stretch", height=300)

    with tabs[2]:
        ctx["targets_editor"]()

    # ---------------------------------------------------------------- extract
    with tabs[3]:
        st.caption("For each ion in each file: the m/z windows are tried in increasing order and the first one "
                   "containing a signal is used. The intensity is the mean over the most stable window of the "
                   "given length. Nothing is excluded for instability; such values are flagged.")
        with st.form("int_settings"):
            c1, c2, c3 = st.columns(3)
            ppm = c1.text_input("m/z windows tried, ± ppm (in this order)", ", ".join(f"{p:g}" for p in ic["ppm_steps"]))
            options = sorted({1, 3, 5, int(ic["min_scans"])})
            min_scans = c1.selectbox("Consecutive scans needed to report a signal", options,
                                     index=options.index(int(ic["min_scans"])))
            length = c2.number_input("Window length, s", value=float(ic["window_length_s"]), min_value=1.0, step=1.0)
            variation = c2.number_input("Flag when windows differ by more than, %",
                                        value=float(ic["variation_flag_pct"]), min_value=1.0, step=5.0,
                                        help="Difference between the highest and lowest window mean of the "
                                             "acquisition, relative to the reported intensity.")
            rsd_flag = c3.number_input("Flag low stability above RSD, %", value=float(ic["rsd_flag_pct"]),
                                       min_value=1.0, step=5.0)
            rel = c3.number_input("Report other peaks from, % of the most intense",
                                  value=float(ic["min_relative_intensity"]) * 100, min_value=1.0, max_value=100.0, step=5.0)
            if st.form_submit_button("Save settings"):
                try:
                    ic.update(ppm_steps=sorted(parse_numbers(ppm)), min_scans=int(min_scans),
                              window_length_s=float(length), variation_flag_pct=float(variation),
                              rsd_flag_pct=float(rsd_flag), min_relative_intensity=float(rel) / 100)
                    save_cfg(root, cfg)
                    st.success("Settings saved.")
                except ValueError as e:
                    st.error(f"Not saved: {e}")
        ready = status["Components"] and len(discover_files(Path(paths["mzml_dir"]))) > 0
        if not ready:
            st.warning("Save the components (tab 1) and provide converted files (tab 2) first.")
        if st.button("Extract intensities", type="primary", disabled=not ready,
                     help="Files already processed with the same settings and targets are reused."):
            bar = st.progress(0.0, text="Starting...")
            buf = io.StringIO()
            try:
                with contextlib.redirect_stdout(buf):
                    cli.cmd_intensity(resolved(root), None,
                                      progress=lambda i, n, f: bar.progress(i / n, text=f"{i} / {n}: {f}"))
                bar.empty()
                st.session_state["int_run_log"] = buf.getvalue()
                st.rerun()
            except (Exception, SystemExit) as e:
                bar.empty()
                st.error(f"Extraction failed: {e}")
                st.session_state["int_run_log"] = buf.getvalue()
        if st.session_state.get("int_run_log"):
            st.code(st.session_state["int_run_log"])

    # ---------------------------------------------------------------- results
    with tabs[4]:
        values = _load(rep / "intensity_values.csv")
        if values.empty:
            st.info("No results yet. Extract the intensities in tab 4.")
            return
        files = _load(rep / "intensity_files.csv")
        peaks = _load(rep / "intensity_peaks.csv")
        matrix = _load(rep / "intensity_matrix.csv")
        values["flags"] = values["flags"].fillna("")
        values["folder"] = values["folder"].fillna("")
        det = values[values["status"] == "detected"]

        m = st.columns(5)
        failed = files[files["processing"] == "failed"]
        m[0].metric("Files", len(files), help=f"{len(failed)} could not be read")
        m[1].metric("Ions detected", det["target_id"].nunique(), help="in at least one file")
        m[2].metric("Intensity values", len(det))
        m[3].metric("Flagged values", int((det["flags"] != "").sum()))
        m[4].metric("Not detected", int((values["status"] == "not detected").sum()),
                    help="ion and file combinations with no signal in any m/z window")
        if len(failed):
            with st.expander(f"{len(failed)} file(s) could not be read"):
                st.dataframe(failed[["file", "error"]], hide_index=True, width="stretch")

        st.subheader("Intensity table")
        fl = st.columns([2, 2, 1])
        folders = sorted(values["folder"].unique())
        pick_folders = fl[0].multiselect("Folders", folders, default=folders,
                                         format_func=lambda x: x or "(top folder)")
        ions = sorted(det["target_id"].unique())
        pick_ions = fl[1].multiselect("Ions", ions, default=ions)
        layout = fl[2].radio("Layout", ["files × ions", "one row per value"])
        if layout == "files × ions":
            keep_cols = ["folder", "name"] + [c for c in pick_ions if c in matrix]
            mx = matrix[matrix["folder"].fillna("").isin(pick_folders)][keep_cols]
            st.dataframe(mx, hide_index=True, width="stretch",
                         column_config={c: st.column_config.NumberColumn(format="%.3e") for c in pick_ions})
            st.caption("Empty cell: ion not detected in that file. Flags are in the 'one row per value' layout.")
            st.download_button("Download this table (.csv)", mx.to_csv(index=False), "intensity_matrix.csv", "text/csv")
        else:
            only_flag = st.checkbox("Only flagged values")
            lv = values[values["folder"].isin(pick_folders) & (values["target_id"].isin(pick_ions)
                                                              | (values["status"] != "detected"))]
            if only_flag:
                lv = lv[lv["flags"] != ""]
            show = ["folder", "name", "target_id", "status", "intensity", "flags", "rsd_pct", "window_variation_pct",
                    "mz_window_ppm", "ppm_error", "mz_observed", "window_start_s", "window_end_s",
                    "scans_in_window", "selected_by"]
            st.dataframe(lv[[c for c in show if c in lv]], hide_index=True, width="stretch",
                         column_config={"intensity": st.column_config.NumberColumn(format="%.3e"),
                                        "rsd_pct": st.column_config.NumberColumn("RSD %", format="%.1f"),
                                        "window_variation_pct": st.column_config.NumberColumn("variation %", format="%.0f"),
                                        "ppm_error": st.column_config.NumberColumn(format="%.2f"),
                                        "mz_observed": st.column_config.NumberColumn(format="%.5f")})
            st.download_button("Download this table (.csv)", lv.to_csv(index=False), "intensity_values.csv", "text/csv")
        with st.expander("Meaning of the flags"):
            st.markdown("\n".join(f"- **{k}**: {v}" for k, v in FLAG_HELP.items() if k in FLAGS))

        # ------------------------------------------------- choice of the peak
        st.subheader("Choice of the peak when several are found")
        st.caption("When the m/z window holds several peaks, the one closest to the theoretical m/z, the most "
                   "intense one and those of at least the set percentage of the most intense are listed. "
                   "Tick the peak to export; unticked groups follow the rule below.")
        rules = {"closest": "Closest to the theoretical m/z", "most_intense": "Most intense"}
        rule = st.radio("Rule for all other ions", list(rules), index=list(rules).index(ic["default_peak"]),
                        format_func=rules.get, horizontal=True)
        multi = peaks[(peaks["status"] == "detected") & (peaks["n_peaks_in_window"] > 1)].copy() if len(peaks) else peaks
        edited = None
        if len(multi):
            multi = multi.sort_values(["file", "target_id", "mz_observed"])
            cols = ["selected", "file", "target_id", "mz_observed", "ppm_error", "intensity",
                    "relative_to_most_intense_pct", "rsd_pct", "is_closest", "is_most_intense"]
            edited = st.data_editor(
                multi[cols].rename(columns={"selected": "use"}), hide_index=True, width="stretch",
                disabled=[c for c in cols if c != "selected"], key="peak_editor",
                column_config={"mz_observed": st.column_config.NumberColumn("m/z", format="%.5f"),
                               "ppm_error": st.column_config.NumberColumn("ppm", format="%.2f"),
                               "intensity": st.column_config.NumberColumn(format="%.3e"),
                               "relative_to_most_intense_pct": st.column_config.NumberColumn("% of most intense", format="%.0f"),
                               "rsd_pct": st.column_config.NumberColumn("RSD %", format="%.1f"),
                               "is_closest": "closest", "is_most_intense": "most intense"})
        else:
            st.info("No ion has more than one peak in its m/z window.")
        if st.button("Apply the choice of peaks"):
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
                save_cfg(root, cfg)
                cli.write_intensity_outputs(resolved(root))
                st.session_state.pop("peak_editor", None)
                st.rerun()
            except ValueError as e:
                st.error(f"Not applied: {e}")

        # ------------------------------------------------------------ detail
        st.subheader("Detail of one ion in one file")
        d1, d2 = st.columns(2)
        det_files = sorted(det["file"].unique())
        if not det_files:
            return
        f_sel = d1.selectbox("File", det_files)
        ion_opts = sorted(det.loc[det["file"] == f_sel, "target_id"])
        t_sel = d2.selectbox("Ion", ion_opts)
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

        tr_path = work / "intensity_traces.parquet"
        g1, g2 = st.columns(2)
        if tr_path.exists():
            tr = pd.read_parquet(tr_path, filters=[("file", "==", f_sel), ("target_id", "==", t_sel)])
            if len(tr):
                order = sorted(tr["mz_observed"].unique())
                tr["peak"] = tr["mz_observed"].map(lambda v: f"m/z {v:.5f}")
                labels = [f"m/z {v:.5f}" for v in order]
                win = pd.DataFrame({"x": [row["window_start_s"]], "x2": [row["window_end_s"]]})
                band = alt.Chart(win).mark_rect(color=BLUE, opacity=0.10).encode(x="x:Q", x2="x2:Q")
                enc = dict(x=alt.X("time_s:Q", title="Acquisition time (s)"),
                           y=alt.Y("intensity:Q", title="Absolute intensity", scale=alt.Scale(zero=True)),
                           tooltip=[alt.Tooltip("peak:N"), alt.Tooltip("time_s:Q", format=".1f"),
                                    alt.Tooltip("intensity:Q", format=".3s")])
                if len(order) > 1:
                    enc["color"] = alt.Color("peak:N", title="Peak", legend=alt.Legend(orient="bottom"),
                                             scale=alt.Scale(domain=labels, range=PEAK_COLORS[:len(labels)]))
                    line = alt.Chart(tr).mark_line(strokeWidth=1.6).encode(**enc)
                else:
                    line = alt.Chart(tr).mark_line(strokeWidth=1.6, color=BLUE).encode(**enc)
                with g1:
                    st.markdown("**Signal over the acquisition**")
                    st.altair_chart((band + line).properties(height=300), width="stretch")
                    st.caption(f"Shaded: the window used ({row['window_start_s']:g}–{row['window_end_s']:g} s), "
                               f"the most stable of its length. Exported peak: m/z {row['mz_observed']:.5f}.")
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
                    lo, hi = w["mean"].min(), w["mean"].max()
                    st.caption(f"One dot per {ic['window_length_s']:g} s window. Lowest window mean {lo:.3e}, "
                               f"highest {hi:.3e}: a spread of {row['window_variation_pct']:.0f} % of the reported "
                               f"intensity (flagged from {ic['variation_flag_pct']:g} %).")

        st.divider()
        if st.button("Prepare archive of all tables"):
            st.session_state["int_zip"] = Path(shutil.make_archive(str(work / "intensity_report"), "zip", rep)).read_bytes()
        if "int_zip" in st.session_state:
            st.download_button("Download all tables (.zip)", st.session_state["int_zip"], "mzstab_intensity.zip")
