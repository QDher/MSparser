"""Folder-based organisation of the intensity extraction.

The user gives a data folder, with or without its subfolders. Every folder
that holds acquisitions is a self-contained unit: an `MSparser` subfolder is
created in it with the converted files, the compound list, the ion list, the
settings used and the results.

    <folder>/
      sample1.RAW ...
      MSparser/
        mzml/            converted files
        components.yaml  compound list of this folder
        config.yaml      settings used (rewritten at each run)
        work/            ion list (targets.csv), cache
        results/         intensity tables of this folder

Tables combining all folders are written to <top folder>/MSparser/results_all.
`<folder>/MSparser/config.yaml` is a normal configuration: the folder can be
processed alone with `msparser intensity -c <folder>/MSparser/config.yaml`.
"""
from __future__ import annotations

import contextlib
import copy
import hashlib
import io
import shutil
from pathlib import Path

import pandas as pd
import yaml

from . import cli
from .config import DEFAULTS, _merge, load_config
from .convert import convert_all, find_raw

SUB = "MSparser"
UNIT_PATHS = {"samples": "samples.csv", "components": "components.yaml", "raw_dir": "..",
              "mzml_dir": "mzml", "work_dir": "work", "report_dir": "results"}
DATA_EXT = (".mzml", ".csv")


# --------------------------------------------------------------- discovery
def _has_data(folder: Path) -> bool:
    try:
        return any(p.suffix.lower() == ".raw" or (p.is_file() and p.suffix.lower() in DATA_EXT)
                   for p in folder.iterdir())
    except OSError:
        return False


def find_units(root: str | Path, include_subfolders: bool = True) -> list[Path]:
    """Folders holding acquisitions (raw, mzML or CSV files): the folder itself
    and, optionally, its subfolders at any depth. The program's own `MSparser`
    folders and Waters `.raw` folders are not entered."""
    root = Path(root)
    if not root.is_dir():
        return []
    found = [root] if _has_data(root) else []
    if include_subfolders:
        def walk(folder: Path):
            try:
                subs = sorted(p for p in folder.iterdir() if p.is_dir())
            except OSError:
                return
            for p in subs:
                if p.name == SUB or p.suffix.lower() == ".raw" or p.name.startswith("."):
                    continue
                if _has_data(p):
                    found.append(p)
                walk(p)
        walk(root)
    return found


def label(folder: Path, root: Path) -> str:
    """A unit as shown to the user: path from the top folder, starting with its name."""
    rel = Path(folder).relative_to(root)
    return str(Path(root.name) / rel) if str(rel) != "." else root.name


def unit_files(folder: Path) -> pd.DataFrame:
    """Converted or directly usable files of one unit: MSparser/mzml first,
    then mzML/CSV files lying in the folder itself."""
    rows, seen = [], set()
    for d in (folder / SUB / "mzml", folder):
        if not d.is_dir():
            continue
        for p in sorted(d.iterdir()):
            if p.is_file() and p.suffix.lower() in DATA_EXT and p.stem not in seen:
                seen.add(p.stem)
                rows.append(dict(file=p.stem, folder="", name=p.stem, path=str(p)))
    return pd.DataFrame(rows, columns=["file", "folder", "name", "path"])


def unit_status(folder: Path) -> dict:
    sub = folder / SUB
    comp = sub / "components.yaml"
    n_comp = len((yaml.safe_load(comp.read_text()) or {}).get("components") or []) if comp.exists() else 0
    values = sub / "results" / "intensity_values.csv"
    return dict(raw_files=len(find_raw(folder)), usable_files=len(unit_files(folder)), compounds=n_comp,
                ions=(sub / "work" / "targets.csv").exists(), results=values.exists())


# ----------------------------------------------------------------- settings
def load_settings(root: Path) -> dict:
    """Settings of a top folder (its MSparser/config.yaml), completed with the defaults."""
    p = Path(root) / SUB / "config.yaml"
    user = (yaml.safe_load(p.read_text()) or {}) if p.exists() else {}
    cfg = _merge(copy.deepcopy(DEFAULTS), user)
    cfg["mode"] = "intensity"
    cfg["paths"] = dict(UNIT_PATHS)
    cfg["intensity"]["folders"] = None
    return cfg


def save_settings(folder: Path, cfg: dict) -> Path:
    """Write the settings into <folder>/MSparser/config.yaml, with the fixed layout of a unit."""
    out = copy.deepcopy(cfg)
    out["mode"] = "intensity"
    out["paths"] = dict(UNIT_PATHS)
    out["intensity"]["folders"] = None
    sub = Path(folder) / SUB
    sub.mkdir(parents=True, exist_ok=True)
    (sub / "config.yaml").write_text(yaml.safe_dump(out, sort_keys=False))
    return sub / "config.yaml"


def unit_cfg(folder: Path, settings: dict, write: bool = False) -> dict:
    """Resolved configuration of one unit. With `write`, the settings are
    saved in the unit first (they are a record of what was used)."""
    sub = Path(folder) / SUB
    if write:
        return load_config(save_settings(folder, settings))
    cfg = copy.deepcopy(settings)
    cfg["paths"] = {k: str((sub / v).resolve()) for k, v in UNIT_PATHS.items()}
    cfg["extraction"]["ppm_windows"] = sorted({float(p) for p in cfg["extraction"]["ppm_windows"]}, reverse=True)
    return cfg


# -------------------------------------------------------------- compound list
def copy_components(src: Path, targets: list[Path]) -> int:
    """Copy the compound list of one folder to other folders (replacing theirs)."""
    source = Path(src) / SUB / "components.yaml"
    n = 0
    for t in targets:
        if Path(t) == Path(src):
            continue
        (Path(t) / SUB).mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, Path(t) / SUB / "components.yaml")
        n += 1
    return n


def _components_hash(path: Path) -> str:
    return hashlib.sha1(path.read_bytes()).hexdigest()


# -------------------------------------------------------------------- running
def run_unit(folder: Path, settings: dict, convert: bool = True, progress=None) -> dict:
    """Convert (if needed) and extract one folder. Never raises: the outcome is
    returned as a row (status, message, counts)."""
    folder = Path(folder)
    row = dict(folder=str(folder), status="", message="", files=0, values=0, flagged=0, not_detected=0, log="")
    buf = io.StringIO()
    try:
        cfg = unit_cfg(folder, settings, write=True)
        sub, work = folder / SUB, Path(cfg["paths"]["work_dir"])
        comp = Path(cfg["paths"]["components"])
        if not comp.exists():
            row.update(status="skipped", message="no compound list in this folder")
            return row
        work.mkdir(parents=True, exist_ok=True)
        with contextlib.redirect_stdout(buf):
            if convert and find_raw(folder):
                try:
                    convert_all(cfg, recursive=False)
                except Exception as e:  # keep going with the files already converted
                    row["message"] = f"conversion: {e}"
            # the ion list follows the compound list: regenerated when the list changed
            mark, digest = work / "components.sha1", _components_hash(comp)
            if not (work / "targets.csv").exists() or not mark.exists() or mark.read_text() != digest:
                cli.cmd_targets(cfg, type("A", (), {"force": True})())
                mark.write_text(digest)
            files = unit_files(folder)
            if files.empty:
                row.update(status="failed", message=row["message"] or "no converted file (mzML or CSV) in this folder")
                return row
            cli.cmd_intensity(cfg, None, progress, files)
        values = pd.read_csv(sub / "results" / "intensity_values.csv")
        det = values[values["status"] == "detected"]
        f = pd.read_csv(sub / "results" / "intensity_files.csv")
        n_failed = int((f["processing"] == "failed").sum())
        row.update(status="done", files=len(f), values=len(det),
                   flagged=int(det["flags"].fillna("").ne("").sum()),
                   not_detected=int((values["status"] == "not detected").sum()))
        if n_failed:
            row["message"] = (row["message"] + "; " if row["message"] else "") + f"{n_failed} file(s) could not be read"
    except Exception as e:
        row.update(status="failed", message=f"{type(e).__name__}: {e}")
    finally:
        row["log"] = buf.getvalue()
    return row


def combine(root: Path, folders: list[Path]) -> dict:
    """Tables over all folders, written to <root>/MSparser/results_all."""
    root = Path(root)
    longs, files = [], []
    for f in folders:
        res = Path(f) / SUB / "results"
        if (res / "intensity_values.csv").exists():
            v = pd.read_csv(res / "intensity_values.csv")
            v.insert(0, "data_folder", label(Path(f), root))
            longs.append(v.drop(columns=[c for c in ("folder",) if c in v]))
        if (res / "intensity_files.csv").exists():
            x = pd.read_csv(res / "intensity_files.csv")
            x.insert(0, "data_folder", label(Path(f), root))
            files.append(x.drop(columns=[c for c in ("folder",) if c in x]))
    out = root / SUB / "results_all"
    if not longs:
        return dict(long=pd.DataFrame(), matrix=pd.DataFrame(), files=pd.DataFrame())
    out.mkdir(parents=True, exist_ok=True)
    long = pd.concat(longs, ignore_index=True)
    allf = pd.concat(files, ignore_index=True) if files else pd.DataFrame()
    det = long[long["status"] == "detected"]
    every = long[["data_folder", "file"]].drop_duplicates()
    if len(allf):
        every = pd.concat([every, allf[["data_folder", "file"]]]).drop_duplicates()

    def wide(col):
        m = det.pivot_table(index=["data_folder", "file"], columns="target_id", values=col, aggfunc="first")
        return every.merge(m.reset_index(), on=["data_folder", "file"], how="left").sort_values(["data_folder", "file"])

    matrix = wide("intensity")
    first = ["data_folder", "file", "target_id", "status", "intensity", "mz_theoretical", "mz_observed", "ppm_error"]
    long = long[[c for c in first if c in long] + [c for c in long.columns if c not in first]]
    long.to_csv(out / "intensity_values_all.csv", index=False)
    matrix.to_csv(out / "intensity_matrix_all.csv", index=False)
    wide("mz_observed").to_csv(out / "intensity_matrix_mz_detected_all.csv", index=False)
    wide("ppm_error").to_csv(out / "intensity_matrix_ppm_error_all.csv", index=False)
    if len(allf):
        allf.to_csv(out / "intensity_files_all.csv", index=False)
    return dict(long=long, matrix=matrix, files=allf)


def run_all(root: Path, folders: list[Path], settings: dict, convert: bool = True, progress=None) -> pd.DataFrame:
    """Process the folders one after another, then write the combined tables.
    `progress(i, n, folder, j, m, file)` reports folder i of n and file j of m."""
    rows = []
    for i, f in enumerate(folders, 1):
        cb = (lambda j, m, name, i=i, f=f: progress(i, len(folders), f, j, m, name)) if progress else None
        if progress:
            progress(i, len(folders), f, 0, 0, "")
        rows.append(run_unit(f, settings, convert, cb))
    combine(root, folders)
    out = pd.DataFrame(rows)
    out.insert(0, "data_folder", [label(Path(f), Path(root)) for f in folders])
    return out
