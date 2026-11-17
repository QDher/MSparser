import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).parent))

from mzstab.cli import main
from mzstab.targets import Component, generate_targets


def _targets():
    comps = [
        Component("caffeine", "C8H10N4O2", (1, 2, -1), "proton"),
        Component("sodium", "Na", (1,), "intrinsic"),
        Component("chloride", "Cl", (-1,), "intrinsic"),
    ]
    return generate_targets(comps).set_index("target_id")


def test_target_masses():
    t = _targets()
    assert t.loc["caffeine [M+H]+", "mz"] == pytest.approx(195.08765, abs=2e-5)
    assert t.loc["caffeine [M+Na]+", "mz"] == pytest.approx(217.06960, abs=2e-5)
    assert t.loc["caffeine [M+2H]2+", "mz"] == pytest.approx(98.04746, abs=2e-5)
    assert t.loc["caffeine [M-H]-", "mz"] == pytest.approx(193.07310, abs=2e-5)
    assert t.loc["caffeine [M+Cl]-", "mz"] == pytest.approx(229.04978, abs=2e-5)
    assert t.loc["sodium [M]+", "mz"] == pytest.approx(22.98922, abs=2e-5)


def test_target_rules():
    t = _targets()
    ions = set(t.index)
    assert {"caffeine [M+H+Na]2+", "caffeine [M+2Na]2+"} <= ions
    assert not any("2M" in i for i in ions)                 # no dimers
    assert "caffeine [M-2H]2-" not in ions                  # -2 not allowed by the user
    assert t.loc["caffeine [M+Na]+", "requires"] == "caffeine;sodium"
    # without sodium entered, no sodium adduct is proposed
    t2 = generate_targets([Component("caffeine", "C8H10N4O2", (1,), "proton")])
    assert list(t2["ion"]) == ["[M+H]+"]
    # 13C isotopologue of caffeine [M+H]+
    assert t.loc["caffeine [M+H]+", "iso_ratio"] == pytest.approx(8 * 0.0107 / 0.9893, rel=1e-6)


@pytest.fixture(scope="module")
def project(tmp_path_factory):
    from synth import make_project

    root = tmp_path_factory.mktemp("proj")
    cfg = make_project(root)
    main(["run", "-c", str(cfg)])
    return root


def test_end_to_end(project):
    rep = project / "report"
    summ = pd.read_csv(rep / "trace_summary.csv")
    s10 = summ[summ["ppm"] == 10].set_index(["file", "target_id"])
    # MS2 scans dropped
    assert pd.read_csv(rep / "scans_dropped.csv")["n_scans"].eq(2).all()
    # Na+ below the scan range; undetected ions flagged
    assert s10.loc[("run02", "sodium [M]+"), "out_of_range"]
    assert s10.loc[("run02", "caffeine [M+2H]2+"), "not_detected"]
    assert s10.loc[("run02", "caffeine [M+H]+"), "analyzed"]
    assert s10.loc[("run02", "caffeine [M+H]+"), "iso_status"] == "ok"
    assert abs(s10.loc[("run02", "caffeine [M+H]+"), "ppm_error_median"]) < 1

    res = pd.read_csv(rep / "trace_results.csv")
    r = res[(res["ppm"] == 10) & (res["axis"] == "clock")].set_index(["file", "target_id"])
    ok = r.loc[("run02", "caffeine [M+H]+")]
    assert not ok["no_stable_window"]
    assert ok["start"] >= 4 and ok["start"] + ok["length"] <= 26.5   # inside the plateau
    assert 3 <= ok["plateau_start_s"] <= 7 and 23 <= ok["plateau_end_s"] <= 27
    assert r.loc[("run08", "caffeine [M+H]+"), "no_stable_window"]     # drifting file
    assert r.loc[("run08", "caffeine [M+H]+"), "trend_rel_slope"] < 0

    gw = pd.read_csv(rep / "global_windows.csv")
    g = gw[(gw["ppm"] == 10) & (gw["axis"] == "clock")].iloc[0]
    assert g["start"] >= 4 and g["start"] + g["length"] <= 26.5
    assert g["pct_all"] >= 60

    fs = pd.read_csv(rep / "file_summary.csv")
    fs = fs[(fs["ppm"] == 10) & (fs["axis"] == "clock")].set_index("file")
    assert (fs["n_ions"] == 4).all() and len(fs) == 8
    assert fs.loc["run02", "n_with_stable_window"] == 4 and fs.loc["run02", "n_stable_in_global_window"] == 4
    assert fs.loc["run08", "n_with_stable_window"] == 0                 # drifting file
    assert "caffeine [M+H]+" in fs.loc["run08", "ions_without_stable_window"]

    grp = pd.read_csv(rep / "group_windows.csv")
    assert (grp["n_acquisitions"] == 8).all()
    assert set(grp["target_id"]) == {"caffeine [M+H]+", "caffeine [M+Na]+", "choline [M]+",
                                     "glucose [M+Na]+"}

    thr = json.loads((rep / "thresholds.json").read_text())
    assert thr["ppm10"]["model"] == "noise" and thr["ppm10"]["a"] > 0

    # interference at +7 ppm: present at 10 ppm, removed at 5 ppm
    cmp_ = pd.read_csv(rep / "window_comparison.csv").set_index(["file", "target_id"])
    assert cmp_.loc[("run01", "caffeine [M+H]+"), "interference_removed_fraction"] > 0.5
    assert cmp_.loc[("run02", "caffeine [M+H]+"), "interference_removed_fraction"] == 0
    assert cmp_.loc[("run02", "caffeine [M+H]+"), "lost_fraction"] < 0.05

    assert len(pd.read_csv(rep / "validation.csv")) > 0
    for f in ("profiles_ppm10.png", "stable_pct_clock_ppm5.png", "rsd_vs_intensity_ppm10.png",
              "mass_error.png", "mixed_model_ppm10.txt", "fixed_window.csv"):
        assert (rep / f).exists(), f


def test_csv_input(project, tmp_path):
    """The same spectra given as CSV (one row per peak) give the same windows."""
    import shutil

    from mzstab.io import read_mzml

    (tmp_path / "mzml").mkdir()
    for f in ("config.yaml", "components.yaml", "samples.csv"):
        shutil.copy(project / f, tmp_path / f)
    for p in sorted((project / "mzml").glob("*.mzML")):
        scans, peaks = read_mzml(p)
        scans = scans[scans["ms_level"] == 1]
        long = peaks.merge(scans[["scan_idx", "time_s"]], on="scan_idx")
        long = long.rename(columns={"time_s": "RT (min)", "mz": "m/z", "intensity": "Intensity"})
        long["RT (min)"] /= 60
        long[["RT (min)", "m/z", "Intensity"]].to_csv(tmp_path / "mzml" / f"{p.stem}.csv", index=False)
    with open(tmp_path / "config.yaml", "a") as fh:
        fh.write("csv: {time_unit: min, polarity: 1}\n")
    main(["run", "-c", str(tmp_path / "config.yaml")])
    cols = ["ppm", "axis", "start", "length", "n_stable", "n_traces"]
    a = pd.read_csv(project / "report" / "global_windows.csv")[cols]
    b = pd.read_csv(tmp_path / "report" / "global_windows.csv")[cols]
    pd.testing.assert_frame_equal(a, b)
    # no scan range in CSV: the out-of-range check is skipped, Na+ is 'not detected'
    s = pd.read_csv(tmp_path / "report" / "trace_summary.csv")
    na = s[s["target_id"] == "sodium [M]+"]
    assert not na["out_of_range"].any() and na["not_detected"].all()


def test_automatic_settings(project, tmp_path):
    """With only the ppm windows configured, the other settings are detected
    from the files and give the same windows as the explicit configuration."""
    import shutil

    shutil.copytree(project / "mzml", tmp_path / "mzml")
    for f in ("components.yaml", "samples.csv"):
        shutil.copy(project / f, tmp_path / f)
    (tmp_path / "config.yaml").write_text(
        "extraction:\n  ppm_windows: [10, 5]\nthresholds: {rsd_factor: auto, slope_limit: auto}\n")
    main(["run", "-c", str(tmp_path / "config.yaml")])
    used = json.loads((tmp_path / "report" / "settings_used.json").read_text())
    assert used["windows"]["lengths_s"] == [3, 5, 8, 10, 15]
    assert used["windows"]["step_s"] == 1.0
    f = used["files"]["run01"]
    assert f["vendor"] == "thermo" and f["polarity"] == 1 and f["n_scans"] == 90
    assert f["scan_range"] == [50.0, 750.0]
    assert f["resolution"] == {"value": 140000.0, "at_mz": 200.0, "scaling": "orbitrap", "source": "file"}
    cols = ["ppm", "axis", "start", "length", "n_stable", "n_traces"]
    pd.testing.assert_frame_equal(pd.read_csv(project / "report" / "global_windows.csv")[cols],
                                  pd.read_csv(tmp_path / "report" / "global_windows.csv")[cols])


def test_manual_windows(project, tmp_path):
    """Manual windows are evaluated like the automatic ones and exported with them."""
    import shutil

    import yaml

    shutil.copytree(project / "work", tmp_path / "work")
    for f in ("components.yaml", "samples.csv"):
        shutil.copy(project / f, tmp_path / f)
    cfg = yaml.safe_load((project / "config.yaml").read_text())
    cfg["manual_windows"] = {"global": [12, 20], "targets": {"caffeine [M+H]+": [0, 6]}}
    (tmp_path / "config.yaml").write_text(yaml.safe_dump(cfg))
    main(["analyze", "-c", str(tmp_path / "config.yaml")])
    c = pd.read_csv(tmp_path / "report" / "window_choice.csv")
    c = c[c["ppm"] == 10].set_index(["target_id", "window"])
    # all ions: automatic windows plus the manual one
    g = c.loc[("choline [M]+", "manual, all ions")]
    assert (g["start"], g["end"], g["chosen_by"]) == (12, 20, "user")
    assert g["pct_stable"] >= 60 and g["n_acquisitions"] == 8
    assert c.loc[("choline [M]+", "automatic, all ions"), "chosen_by"] == "automatic"
    assert ("choline [M]+", "automatic, this ion") in c.index
    # the weakest ion is too noisy for its slope to be assessed in 8 s
    assert c.loc[("glucose [M+Na]+", "manual, all ions"), "pct_stable"] < g["pct_stable"]
    # the window set for one ion takes precedence; 0-6 s is the rise: not stable
    m = c.loc[("caffeine [M+H]+", "manual, this ion")]
    assert (m["start"], m["end"]) == (0, 6) and m["pct_stable"] == 0
    assert ("caffeine [M+H]+", "manual, all ions") not in c.index
    t = pd.read_csv(tmp_path / "report" / "manual_window_traces.csv")
    assert len(t[t["ppm"] == 10]) == 32 and {"mean", "rsd", "rel_slope", "stable"} <= set(t.columns)
    # removing the manual windows removes their export
    cfg["manual_windows"] = {"global": None, "targets": {}}
    (tmp_path / "config.yaml").write_text(yaml.safe_dump(cfg))
    main(["analyze", "-c", str(tmp_path / "config.yaml")])
    assert not (tmp_path / "report" / "manual_window_traces.csv").exists()
    assert set(pd.read_csv(tmp_path / "report" / "window_choice.csv")["chosen_by"]) == {"automatic"}


def test_pair_targets():
    comps = [
        Component("caffeine", "C8H10N4O2", (1,), "proton"),
        Component("glucose", "C6H12O6", (1, -1), "proton"),
        Component("calcium", "Ca", (2,), "intrinsic"),
        Component("chloride", "Cl", (-1,), "intrinsic"),
    ]
    t = generate_targets(comps).set_index("target_id")
    p = t.loc["caffeine + glucose [M1+M2+H]+"]
    assert p["mz"] == pytest.approx(194.080376 + 180.063388 + 1.007276, abs=2e-5)
    assert p["type"] == "pair" and p["requires"] == "caffeine;glucose"
    assert "caffeine + glucose [M1+M2-H]-" in t.index          # glucose allows -1
    assert "caffeine + glucose [M1+M2+2H]2+" in t.index        # both allow +1
    assert "caffeine + glucose [M1+M2-2H]2-" not in t.index    # caffeine does not allow -1
    assert t.loc["calcium + chloride [M1+M2]+", "mz"] == pytest.approx(39.962591 + 34.968853 - 0.000549, abs=2e-5)
    assert t.loc["glucose [M+Cl]-", "type"] == "adduct" and t.loc["caffeine [M+H]+", "type"] == "ion"
    assert not any("caffeine + caffeine" in i for i in t.index)  # no component with itself
    assert (generate_targets(comps, pairs=False)["type"] != "pair").all()


def test_convert_reports_and_fails_clearly(tmp_path, capsys):
    from mzstab.config import load_config
    from mzstab.convert import convert_all

    (tmp_path / "config.yaml").write_text("vendor: thermo\n")
    cfg = load_config(tmp_path / "config.yaml")
    with pytest.raises(FileNotFoundError, match="raw folder not found"):
        convert_all(cfg)
    (tmp_path / "raw").mkdir()
    with pytest.raises(FileNotFoundError, match="no .raw file"):
        convert_all(cfg)
    (tmp_path / "raw" / "runA.RAW").write_text("x")       # upper-case extension
    with pytest.raises(FileNotFoundError, match="converter not found"):
        convert_all(cfg)

    def fake(name, body):
        """A stand-in converter written in Python, so that it runs on every system."""
        p = tmp_path / f"{name}.py"
        p.write_text("import sys, pathlib\n" + body)
        return f'"{sys.executable}" "{p}"'

    # a converter that exits normally without writing the mzML is an error
    cfg["convert"]["thermo_parser"] = fake("silent", "print('license missing')\n")
    with pytest.raises(RuntimeError, match="did not produce runA.mzML(.|\n)*license missing"):
        convert_all(cfg)
    # a working converter: -i=<raw> -o=<dir>
    cfg["convert"]["thermo_parser"] = fake(
        "good",
        "i = pathlib.Path(sys.argv[1][3:]); o = pathlib.Path(sys.argv[2][3:])\n"
        "(o / (i.stem + '.mzML')).write_text('data')\n")
    convert_all(cfg)
    assert (tmp_path / "mzml" / "runA.mzML").exists()
    assert "runA.RAW -> runA.mzML" in capsys.readouterr().out
    convert_all(cfg)
    assert "already converted" in capsys.readouterr().out


def test_converter_setting_may_be_a_folder(tmp_path):
    from mzstab.convert import _program

    sub = tmp_path / "ThermoRawFileParser-win" / "inner"
    sub.mkdir(parents=True)
    with pytest.raises(FileNotFoundError, match="is a folder and contains no ThermoRawFileParser.exe"):
        _program(str(tmp_path), "thermo_parser")
    exe = sub / "ThermoRawFileParser.exe"
    exe.write_text("x")
    assert _program(str(tmp_path), "thermo_parser") == [str(exe)]      # found in a subfolder
    assert _program(f'"{exe}"', "thermo_parser") == [str(exe)]          # quoted full path


def test_no_stable_window_and_few_traces(project, tmp_path):
    """The analysis completes when no window is stable, and with one or two traces."""
    import shutil

    import yaml

    # (a) thresholds that no window can meet
    a = tmp_path / "a"
    shutil.copytree(project / "work", a / "work")
    for f in ("components.yaml", "samples.csv"):
        shutil.copy(project / f, a / f)
    cfg = yaml.safe_load((project / "config.yaml").read_text())
    cfg["thresholds"] = {"slope_limit": 1e-9}
    (a / "config.yaml").write_text(yaml.safe_dump(cfg))
    main(["analyze", "-c", str(a / "config.yaml")])
    r = pd.read_csv(a / "report" / "trace_results.csv")
    assert r["no_stable_window"].all()
    assert not (a / "report" / "global_windows.csv").exists()

    # (b) two acquisitions, one component: two traces per ppm setting
    b = tmp_path / "b"
    (b / "mzml").mkdir(parents=True)
    for f in ("run01", "run02"):
        shutil.copy(project / "mzml" / f"{f}.mzML", b / "mzml")
    (b / "samples.csv").write_text("file\nrun01\nrun02\n")
    (b / "components.yaml").write_text(
        "components:\n  - {name: caffeine, formula: C8H10N4O2, charges: [1], mode: molecule}\n")
    (b / "config.yaml").write_text("extraction:\n  ppm_windows: [10, 5]\n")
    main(["run", "-c", str(b / "config.yaml")])
    r = pd.read_csv(b / "report" / "trace_results.csv")
    assert len(r[(r["ppm"] == 10) & (r["axis"] == "clock")]) == 2
    thr = json.loads((b / "report" / "thresholds.json").read_text())["ppm10"]
    assert thr["n_traces"] == 2 and thr["model"] == "constant" and thr["slope_limit_derived"] > 0
    # by default the limits are derived from the data
    assert thr["rsd_factor_source"] == thr["slope_limit_source"] == "data"
    assert thr["rsd_factor"] == thr["rsd_factor_derived"] and thr["slope_limit"] == thr["slope_limit_derived"]
    assert not r.loc[(r["ppm"] == 10) & (r["axis"] == "clock"), "no_stable_window"].all()


def test_app_starts(tmp_path):
    """The interface loads on an empty project without raising."""
    AppTest = pytest.importorskip("streamlit.testing.v1").AppTest
    app = Path(__file__).parent.parent / "mzstab" / "app.py"
    old = sys.argv
    sys.argv = ["app.py", str(tmp_path)]
    try:
        at = AppTest.from_file(str(app), default_timeout=120).run()
    finally:
        sys.argv = old
    assert not at.exception, [e.value for e in at.exception]
    assert any(t.label == "Folder containing the raw files" for t in at.text_input)
