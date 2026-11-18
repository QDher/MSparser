import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).parent))

from mzstab.cli import main
from mzstab.intensity import group_peaks

CAF = "caffeine [M+H]+"


def test_group_peaks():
    rng = np.random.default_rng(0)
    c = 195.0877
    one = c * (1 + rng.normal(0, 1e-6, 80))
    assert len(set(group_peaks(one, c, 1.0))) == 1
    two = np.r_[one, c * (1 + 7e-6 + rng.normal(0, 1e-6, 80))]
    lab = group_peaks(two, c, 1.0)
    assert len(set(lab)) == 2 and len(set(lab[:80])) == 1 and len(set(lab[80:])) == 1


@pytest.fixture(scope="module")
def batch(tmp_path_factory):
    """Files organised in subfolders, with different situations per file."""
    from synth import COMPONENTS, write_file

    root = tmp_path_factory.mktemp("batch")
    for sub in ("mzml/day1", "mzml/day2/series_b"):
        (root / sub).mkdir(parents=True)
    write_file(root / "mzml/day1/a_plain.mzML", seed=1)
    write_file(root / "mzml/day1/b_offset8ppm.mzML", seed=2, mz_offset_ppm=8)
    write_file(root / "mzml/day2/c_drift.mzML", seed=3, drift=-4)
    write_file(root / "mzml/day2/series_b/d_small_interf.mzML", seed=4, interference=True)
    write_file(root / "mzml/day2/series_b/e_big_interf.mzML", seed=5, interference=True,
               interference_intensity=3e7)
    (root / "mzml/day2/broken.mzML").write_text("not an mzML file")
    (root / "components.yaml").write_text(COMPONENTS)
    (root / "config.yaml").write_text("mode: intensity\n")
    main(["intensity", "-c", str(root / "config.yaml")])
    return root


def test_batch_outputs(batch, capsys):
    rep = batch / "report"
    files = pd.read_csv(rep / "intensity_files.csv").set_index("file")
    assert len(files) == 6 and files.loc["day2/broken", "processing"] == "failed"
    assert files.loc["day2/series_b/d_small_interf", "folder"] == "day2/series_b"
    assert (files.drop("day2/broken")["processing"] == "processed").all()

    v = pd.read_csv(rep / "intensity_values.csv", keep_default_na=False, na_values=[""]).set_index(["file", "target_id"])
    a = v.loc[("day1/a_plain", CAF)]
    assert a["status"] == "detected" and a["mz_window_ppm"] == 5
    assert a["intensity"] == pytest.approx(1e7, rel=0.03)            # plateau value, not the rise or decay
    assert 5 <= a["window_start_s"] and a["window_end_s"] <= 25.5
    assert a["window_end_s"] - a["window_start_s"] == pytest.approx(10)
    assert a["rsd_pct"] < 4 and pd.isna(a["flags"])
    # rise and decay are inside the acquisition: the 10 s windows differ, which is reported
    assert v.loc[("day1/a_plain", CAF), "n_windows"] > 10

    # adaptive m/z window: 8 ppm off -> not found at 5 ppm, found at 10 ppm
    b = v.loc[("day1/b_offset8ppm", CAF)]
    assert b["mz_window_ppm"] == 10 and b["ppm_error"] == pytest.approx(8, abs=0.5)
    assert "wide_mz_window" in b["flags"] and b["intensity"] == pytest.approx(1e7, rel=0.03)

    # drifting signal: still reported, flagged for variation between windows
    c = v.loc[("day2/c_drift", CAF)]
    assert c["status"] == "detected" and "window_variation" in c["flags"] and c["window_variation_pct"] >= 25

    # ions absent from the data are listed as not detected, with no intensity
    nd = v.loc[("day1/a_plain", "caffeine [M+2H]2+")]
    assert nd["status"] == "not detected" and pd.isna(nd["intensity"])
    assert v.loc[("day1/a_plain", "sodium [M]+"), "status"] == "out of scan range"

    m = pd.read_csv(rep / "intensity_matrix.csv").set_index("file")
    assert m.loc["day1/a_plain", CAF] == pytest.approx(a["intensity"])
    assert pd.isna(m.loc["day2/broken", CAF]) and len(m) == 6
    assert (rep / "intensity_window_means.csv").exists()


def test_candidate_peaks_and_selection(batch):
    """A second peak in the m/z window is reported if it is the most intense or
    at least 10 % of it; the user's choice decides which one is exported."""
    import yaml

    cfg = {"mode": "intensity", "intensity": {"ppm_steps": [10, 20]}}
    (batch / "config.yaml").write_text(yaml.safe_dump(cfg))
    main(["intensity", "-c", str(batch / "config.yaml")])
    peaks = pd.read_csv(batch / "report" / "intensity_peaks.csv")
    small = peaks[(peaks["file"] == "day2/series_b/d_small_interf") & (peaks["target_id"] == CAF)]
    assert len(small) == 1                                   # 3 % of the main peak: not reported
    big = peaks[(peaks["file"] == "day2/series_b/e_big_interf") & (peaks["target_id"] == CAF)].sort_values("mz_observed")
    assert len(big) == 2 and big["multiple_peaks"].all()
    main_pk, interf = big.iloc[0], big.iloc[1]
    assert main_pk["is_closest"] and not main_pk["is_most_intense"] and main_pk["selected"]
    assert interf["is_most_intense"] and interf["ppm_error"] == pytest.approx(7, abs=0.5)
    assert interf["intensity"] == pytest.approx(3e7, rel=0.05)
    assert main_pk["relative_to_most_intense_pct"] == pytest.approx(33, abs=4)

    # the user picks the other peak for that file
    from mzstab.cli import write_intensity_outputs
    from mzstab.config import load_config

    pd.DataFrame([dict(file=interf["file"], target_id=CAF, mz_observed=interf["mz_observed"])]).to_csv(
        batch / "work" / "peak_selection.csv", index=False)
    out = write_intensity_outputs(load_config(batch / "config.yaml"))
    row = out["long"].set_index(["file", "target_id"]).loc[(interf["file"], CAF)]
    assert row["intensity"] == pytest.approx(3e7, rel=0.05) and row["selected_by"] == "user"
    other = out["long"].set_index(["file", "target_id"]).loc[("day1/a_plain", CAF)]
    assert other["selected_by"] == "rule: closest m/z"

    # default rule "most intense" for everything else
    (batch / "work" / "peak_selection.csv").unlink()
    cfg["intensity"]["default_peak"] = "most_intense"
    (batch / "config.yaml").write_text(yaml.safe_dump(cfg))
    out = write_intensity_outputs(load_config(batch / "config.yaml"))
    assert out["long"].set_index(["file", "target_id"]).loc[(interf["file"], CAF), "intensity"] == pytest.approx(3e7, rel=0.05)


def test_cache_reuses_unchanged_files(batch, capsys):
    main(["intensity", "-c", str(batch / "config.yaml")])
    out = capsys.readouterr().out
    assert "0 processed, 5 unchanged (reused), 1 failed" in out


def test_interface_of_intensity_mode(batch):
    AppTest = pytest.importorskip("streamlit.testing.v1").AppTest
    app = Path(__file__).parent.parent / "mzstab" / "app.py"
    old = sys.argv
    sys.argv = ["app.py", str(batch)]
    try:
        at = AppTest.from_file(str(app), default_timeout=180).run()
        assert not at.exception, [e.value for e in at.exception]
        assert [t.label for t in at.tabs] == ["1. Components", "2. Files", "3. Targets", "4. Extract", "5. Results"]
        at = [b for b in at.button if b.label == "Extract intensities"][0].click().run()
        assert not at.exception and not at.error, [e.value for e in at.error]
        assert {m.label: m.value for m in at.metric}["Files"] == "6"
        at = [b for b in at.button if b.label == "Apply the choice of peaks"][0].click().run()
        assert not at.exception and not at.error
    finally:
        sys.argv = old


def test_recursive_conversion(tmp_path):
    """Raw files in subfolders are converted into the same folder tree."""
    from mzstab.config import load_config
    from mzstab.convert import convert_all, find_raw

    (tmp_path / "raw" / "day1" / "set_a").mkdir(parents=True)
    (tmp_path / "raw" / "top.raw").write_text("x")
    (tmp_path / "raw" / "day1" / "set_a" / "deep.RAW").write_text("x")
    (tmp_path / "raw" / "day1" / "waters.raw").mkdir()                # a Waters acquisition is a folder
    (tmp_path / "raw" / "day1" / "waters.raw" / "_FUNC001.DAT").write_text("x")
    assert [p.name for p in find_raw(tmp_path / "raw", recursive=True)] == ["deep.RAW", "waters.raw", "top.raw"]
    conv = tmp_path / "conv.py"
    conv.write_text("import sys, pathlib\ni = pathlib.Path(sys.argv[1][3:]); o = pathlib.Path(sys.argv[2][3:])\n"
                    "(o / (i.stem + '.mzML')).write_text('data')\n")
    (tmp_path / "config.yaml").write_text("vendor: thermo\nmode: intensity\n")
    cfg = load_config(tmp_path / "config.yaml")
    cfg["convert"]["thermo_parser"] = f'"{sys.executable}" "{conv}"'
    convert_all(cfg, recursive=True)
    assert (tmp_path / "mzml" / "top.mzML").exists()
    assert (tmp_path / "mzml" / "day1" / "set_a" / "deep.mzML").exists()
