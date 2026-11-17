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
    (tmp_path / "config.yaml").write_text("extraction:\n  ppm_windows: [10, 5]\n")
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
