"""Vendor raw -> centroided mzML. Thin wrappers around external converters.

ThermoRawFileParser centroids by default. msconvert needs the vendor libraries
(Windows, or the ProteoWizard Docker image under Wine).
"""
from __future__ import annotations

import os
import shlex
import shutil
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


PROGRAM_NAMES = {
    "thermo_parser": ("ThermoRawFileParser.exe", "ThermoRawFileParser"),
    "msconvert": ("msconvert.exe", "msconvert"),
}


def _program(spec: str, key: str = "thermo_parser") -> list[str]:
    """Converter setting -> command prefix. Accepts a program name, the full
    path of the program (spaces allowed), the folder that contains it, or a
    prefixed call such as 'mono /x/ThermoRawFileParser.exe'."""
    spec = str(spec).strip()
    whole = spec[1:-1] if spec.count('"') == 2 and spec.startswith('"') and spec.endswith('"') else spec
    if Path(whole).is_dir():
        # a folder was given: look for the program inside it (and its subfolders)
        folder = Path(whole)
        for name in PROGRAM_NAMES[key]:
            hits = sorted(p for p in folder.rglob(name) if p.is_file())
            if hits:
                return [str(hits[0])]
        raise FileNotFoundError(
            f"'{folder}' is a folder and contains no {PROGRAM_NAMES[key][0]}. "
            f"Give the full path of the program file, for example {folder / PROGRAM_NAMES[key][0]}"
        )
    if Path(whole).exists() or shutil.which(whole):
        return [whole]
    parts = shlex.split(spec, posix=os.name != "nt")
    return [x.strip('"') for x in parts] or [spec]


def find_raw(raw_dir: Path, recursive: bool = False) -> list[Path]:
    """Raw acquisitions in a folder (.raw in any letter case; Waters: folders).
    With `recursive`, subfolders are searched too."""
    if not raw_dir.exists():
        raise FileNotFoundError(f"raw folder not found: {raw_dir}")
    if not recursive:
        return sorted(p for p in raw_dir.iterdir() if p.suffix.lower() == ".raw")
    found = []

    def walk(folder: Path):
        for p in sorted(folder.iterdir()):
            if p.suffix.lower() == ".raw":
                found.append(p)          # a Waters .raw folder is one acquisition: do not enter it
            elif p.is_dir():
                walk(p)

    walk(raw_dir)
    return found


def convert_all(cfg: dict, dry_run: bool = False, recursive: bool = False, folders=None) -> list[list[str]]:
    """Convert every raw acquisition without an mzML. Prints one line per file
    and raises if nothing can be converted or a conversion gives no mzML.
    With `recursive`, subfolders are converted too and the folder tree is
    reproduced in the output folder."""
    raw_dir = Path(cfg["paths"]["raw_dir"])
    out_dir = Path(cfg["paths"]["mzml_dir"])
    raws = find_raw(raw_dir, recursive)
    if folders is not None:  # only the chosen folders (relative to the raw folder)
        keep = set(folders)
        raws = [r for r in raws if ("" if r.parent == raw_dir else r.parent.relative_to(raw_dir).as_posix()) in keep]
    if not raws:
        raise FileNotFoundError(
            f"no .raw file in {raw_dir}. Copy the raw files there, or set the raw folder in the settings."
        )
    if not dry_run:
        out_dir.mkdir(parents=True, exist_ok=True)
    cmds, done, skipped = [], 0, 0
    for raw in raws:
        sub = raw.parent.relative_to(raw_dir)
        dest = out_dir / sub
        target = dest / f"{raw.stem}.mzML"
        shown = (sub / raw.name).as_posix()
        if target.exists():
            skipped += 1
            print(f"{shown}: already converted")
            continue
        # Waters raw data are folders, Thermo raw data are files
        vendor = cfg["vendor"] if cfg["vendor"] != "auto" else ("waters" if raw.is_dir() else "thermo")
        cmd = (thermo_command if vendor == "thermo" else waters_command)(cfg, raw, dest)
        key = "thermo_parser" if vendor == "thermo" else "msconvert"
        cmd = _program(cfg["convert"][key], key) + cmd[1:]
        cmds.append(cmd)
        if dry_run:
            continue
        dest.mkdir(parents=True, exist_ok=True)
        if not (Path(cmd[0]).exists() or shutil.which(cmd[0])):
            name = "ThermoRawFileParser" if vendor == "thermo" else "msconvert (ProteoWizard)"
            raise FileNotFoundError(
                f"converter not found: '{cmd[0]}'. Install {name} and give the full path of its program "
                f"in Settings > Advanced settings > Converters."
            )
        try:
            r = subprocess.run(cmd, capture_output=True, text=True)
        except OSError as e:
            raise RuntimeError(
                f"the converter could not be started: {cmd[0]}\n{e}\n"
                "Check that the converter setting is the program file itself (ending in .exe on Windows), "
                "not a folder or another file. If the program was downloaded, Windows may block it: "
                "right-click the .exe > Properties > tick 'Unblock', or allow it in the antivirus."
            ) from e
        if r.returncode != 0 or not target.exists():
            tail = "\n".join((r.stdout + "\n" + r.stderr).strip().splitlines()[-15:])
            raise RuntimeError(
                f"{shown}: the converter did not produce {target.name} "
                f"(exit code {r.returncode}).\nCommand: {' '.join(cmd)}\nConverter output:\n{tail}"
            )
        done += 1
        print(f"{shown} -> {target.name}")
    if not dry_run:
        print(f"{done} file(s) converted, {skipped} already converted, in {out_dir}")
    return cmds
