"""Vendor raw -> centroided mzML. Thin wrappers around external converters.

ThermoRawFileParser centroids by default. msconvert needs the vendor libraries
(Windows, or the ProteoWizard Docker image under Wine).
"""
from __future__ import annotations

import subprocess
from pathlib import Path


def thermo_command(cfg: dict, raw: Path, out_dir: Path) -> list[str]:
    return [cfg["convert"]["thermo_parser"], f"-i={raw}", f"-o={out_dir}", "-f=2"]


def waters_command(cfg: dict, raw: Path, out_dir: Path) -> list[str]:
    c = cfg["convert"]
    cmd = [c["msconvert"], str(raw), "--mzML", "-o", str(out_dir)]
    if c["lockmass_mz"] is not None:
        f = f"lockmassRefiner mz={c['lockmass_mz']}"
        if c["lockmass_mz_neg"] is not None:
            f += f" mzNegIons={c['lockmass_mz_neg']}"
        f += f" tol={c['lockmass_tol_da']}"
        cmd += ["--filter", f]
    cmd += ["--filter", "peakPicking vendor msLevel=1-"]
    return cmd


def convert_all(cfg: dict, dry_run: bool = False) -> list[list[str]]:
    raw_dir = Path(cfg["paths"]["raw_dir"])
    out_dir = Path(cfg["paths"]["mzml_dir"])
    if not dry_run:
        out_dir.mkdir(parents=True, exist_ok=True)
    pattern = "*.raw"  # Thermo: files; Waters: directories
    cmds = []
    for raw in sorted(raw_dir.glob(pattern)):
        if (out_dir / f"{raw.stem}.mzML").exists():
            continue
        # Waters raw data are folders, Thermo raw data are files
        vendor = cfg["vendor"] if cfg["vendor"] != "auto" else ("waters" if raw.is_dir() else "thermo")
        cmd = (thermo_command if vendor == "thermo" else waters_command)(cfg, raw, out_dir)
        cmds.append(cmd)
        if not dry_run:
            subprocess.run(cmd, check=True)
    return cmds
