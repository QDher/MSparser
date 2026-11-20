"""Interface parts used by both modes: compound list and ion list editors."""
from __future__ import annotations

import argparse
import contextlib
import io
from pathlib import Path

import pandas as pd
import streamlit as st
import yaml

from msparser import cli
from msparser.convert import install_dir, locate_converter
from msparser.targets import Component, generate_targets, normalize_mode

CHARGES = ["+1", "+2", "+3", "-1", "-2", "-3"]
# labels shown in the interface for the two charge modes stored in components.yaml
MODE_LABELS = {"proton": "Molecule (gains or loses H⁺)", "intrinsic": "Ion (already charged)"}
MODE_VALUES = {v: k for k, v in MODE_LABELS.items()}


def read_components(path: Path) -> list[dict]:
    if not Path(path).exists():
        return []
    return (yaml.safe_load(Path(path).read_text()) or {}).get("components", []) or []


def components_editor(path: Path, pairs: bool, key: str, next_step: str = "the Targets tab") -> bool:
    """Editable compound list stored in `path`. Returns True when saved."""
    st.caption("Enter every injected molecule or ion as a neutral formula (for example Na, not Na⁺), choose the "
               "charge states that are possible in the 'Charges' list, and state what the component is:")
    st.markdown(
        f"- **{MODE_LABELS['proton']}**: a neutral molecule that is protonated or deprotonated in the source, "
        "e.g. [M+H]⁺, [M−H]⁻.\n"
        f"- **{MODE_LABELS['intrinsic']}**: a species that is an ion by itself (Na, K, NH4, Cl, a quaternary "
        "ammonium…). No H⁺ is added or removed. **Choose this for an ion to obtain its adducts with the "
        "molecules**: sodium with +1 gives [M+Na]⁺ and [M−2H+Na]⁻; chloride with −1 gives [M+Cl]⁻.\n"
        "- Ions made of two different components are also proposed: two molecules sharing H⁺ "
        "([M1+M2+H]⁺), or a cation with an anion.")
    rows = [{"name": d["name"], "formula": d["formula"],
             "charges": [f"{int(c):+d}" for c in d.get("charges", [])],
             "mode": MODE_LABELS[normalize_mode(d.get("mode", "proton"))]} for d in read_components(path)]
    comp_df = pd.DataFrame(rows, columns=["name", "formula", "charges", "mode"])
    edited = st.data_editor(
        comp_df, num_rows="dynamic", width="stretch", key=f"{key}_editor",
        column_config={
            "name": st.column_config.TextColumn("Name", required=True),
            "formula": st.column_config.TextColumn("Neutral formula", required=True),
            "charges": st.column_config.MultiselectColumn(
                "Charges", options=CHARGES, default=[], width="medium",
                help="Pick every charge state the component can take, from −3 to +3"),
            "mode": st.column_config.SelectboxColumn("The component is a…", options=list(MODE_VALUES),
                                                     default=MODE_LABELS["proton"], required=True, width="medium"),
        })
    if st.button("Save components", type="primary", key=f"{key}_save"):
        try:
            out, comps = [], []
            for r in edited.dropna(subset=["name", "formula"]).to_dict("records"):
                picked = r.get("charges")
                picked = [] if picked is None else list(picked)
                charges = sorted({int(str(c)) for c in picked}, key=lambda v: (v < 0, abs(v)))
                name, mode = str(r["name"]).strip(), MODE_VALUES.get(r.get("mode"), "proton")
                if not charges:
                    raise ValueError(f"{name}: choose at least one charge")
                c = Component(name, str(r["formula"]).strip(), tuple(charges), mode)
                c.composition()  # validates the formula
                comps.append(c)
                out.append({"name": name, "formula": c.formula, "charges": charges, "mode": mode})
            if len({c.name for c in comps}) != len(comps):
                raise ValueError("component names must be unique")
            if comps:
                generate_targets(comps, pairs)
            Path(path).parent.mkdir(parents=True, exist_ok=True)
            Path(path).write_text(yaml.safe_dump({"components": out}, sort_keys=False))
            st.success(f"{len(out)} components saved. Generate the ions again in {next_step}.")
            return True
        except Exception as e:
            st.error(f"Not saved: {e}")
    return False


def targets_editor(cfg: dict, key: str) -> None:
    """Proposed ions of one project or folder (`cfg` holds its resolved paths)."""
    work, comp_path = Path(cfg["paths"]["work_dir"]), Path(cfg["paths"]["components"])
    st.caption("Ions proposed from the components. Untick 'use' for ions that should not be followed.")
    tp = work / "targets.csv"
    if st.button("Generate the ions from the components" if not tp.exists()
                 else "Regenerate the ions (discards the current selection)", key=f"{key}_gen",
                 disabled=not comp_path.exists()):
        buf = io.StringIO()
        try:
            work.mkdir(parents=True, exist_ok=True)
            with contextlib.redirect_stdout(buf):
                cli.cmd_targets(cfg, argparse.Namespace(force=True))
        except (Exception, SystemExit) as e:
            st.error(f"Ions not generated: {e}")
    if not comp_path.exists():
        st.info("Save a compound list first.")
    if not tp.exists():
        return
    t = pd.read_csv(tp)
    t["use"] = t["use"].astype(bool)
    if "type" in t:
        n = t["type"].value_counts()
        st.caption(f"{len(t)} ions: {n.get('ion', 0)} of single components, {n.get('adduct', 0)} adducts of a "
                   f"molecule with entered ions, {n.get('pair', 0)} made of two components. Positive ions "
                   "are followed in positive-mode acquisitions, negative ions in negative-mode ones.")
        cps = read_components(comp_path)
        if cps and not any(normalize_mode(c.get("mode", "proton")) == "intrinsic" for c in cps):
            st.warning(f"No component is set as '{MODE_LABELS['intrinsic']}', so no adduct with an entered ion "
                       "(such as [M+Na]⁺) is proposed. If a component is an ion (Na, K, Cl…), set its "
                       "type to this in the compound list and regenerate the ions.")
    t_edit = st.data_editor(
        t, width="stretch", hide_index=True, key=f"{key}_editor",
        disabled=[c for c in t.columns if c != "use"],
        column_order=["use", "type", "component", "ion", "charge", "mz", "iso_mz", "iso_ratio", "requires"],
        column_config={"mz": st.column_config.NumberColumn("m/z", format="%.5f"),
                       "iso_mz": st.column_config.NumberColumn("isotope m/z", format="%.5f"),
                       "iso_ratio": st.column_config.NumberColumn("isotope ratio", format="%.4f"),
                       "requires": "requires components"})
    if st.button("Save the selection of ions", type="primary", key=f"{key}_save"):
        t_edit.assign(use=t_edit["use"].astype(int)).to_csv(tp, index=False)
        st.success(f"{int(t_edit['use'].sum())} of {len(t_edit)} ions selected.")


def converter_panel(cfg: dict, save, prefix: str) -> None:
    """Where the raw-file converter is, with its status. `cfg` is the settings
    dictionary holding `convert`; `save()` stores it after a change."""
    thermo, msc = locate_converter(cfg, "thermo_parser"), locate_converter(cfg, "msconvert")
    if thermo or msc:
        label = "Converter found: " + ", ".join(n for n, f in (("ThermoRawFileParser", thermo), ("msconvert", msc)) if f)
    else:
        label = "Converter not found — open to set it"
    with st.expander(label, expanded=not (thermo or msc)):
        st.markdown(
            "Raw files are converted by a separate free program: **ThermoRawFileParser** (Thermo) or "
            "**msconvert** from ProteoWizard (Waters). It is found automatically in either case:\n"
            f"- its folder is placed inside the folder of this program: `{install_dir()}` "
            "(for example unzip ThermoRawFileParser there);\n"
            "- or its location is given below, as the program file or the folder that contains it.")
        st.write(("✅ ThermoRawFileParser: " + thermo) if thermo else "⬜ ThermoRawFileParser: not found")
        st.write(("✅ msconvert: " + msc) if msc else "⬜ msconvert: not found (only needed for Waters files)")
        c1, c2 = st.columns(2)
        t_in = c1.text_input("ThermoRawFileParser: program or folder", cfg["convert"]["thermo_parser"],
                             key=f"{prefix}_thermo", help="Example: C:\\Tools\\ThermoRawFileParser. "
                             "Leave 'ThermoRawFileParser' to search the folder of this program.")
        m_in = c2.text_input("msconvert: program or folder (Waters)", cfg["convert"]["msconvert"],
                             key=f"{prefix}_msconvert", help="Example: C:\\Program Files\\ProteoWizard")

        def _clean(v: str) -> str:
            v = v.strip()  # a single quoted path loses its quotes; a quoted multi-part command is kept
            return v[1:-1] if v.count('"') == 2 and v.startswith('"') and v.endswith('"') else v

        new = {"thermo_parser": _clean(t_in) or "ThermoRawFileParser", "msconvert": _clean(m_in) or "msconvert"}
        if any(new[k] != cfg["convert"][k] for k in new):
            cfg["convert"].update(new)
            save()
            st.rerun()
