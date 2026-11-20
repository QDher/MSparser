"""Folder-based intensity extraction: each data folder keeps its own files."""
import sys
from pathlib import Path

import pandas as pd
import pytest
import yaml

sys.path.insert(0, str(Path(__file__).parent))

from msparser import folders as fo
from msparser.cli import main

CAF = "caffeine [M+H]+"
ONLY_CAFFEINE = "components:\n  - {name: caffeine, formula: C8H10N4O2, charges: [1], mode: molecule}\n"


@pytest.fixture(scope="module")
def tree(tmp_path_factory):
    """root/a.RAW, root/sub1/b.RAW (raw files with a stand-in converter),
    root/sub1/deep/c.mzML (already converted), root/empty/."""
    from synth import COMPONENTS, write_file

    base = tmp_path_factory.mktemp("data")
    root, bank = base / "root", base / "bank"
    for d in (root / "sub1" / "deep", root / "empty", bank):
        d.mkdir(parents=True)
    write_file(bank / "a.mzML", seed=1)
    write_file(bank / "b.mzML", seed=2, mz_offset_ppm=8)
    (root / "a.RAW").write_text("raw")
    (root / "sub1" / "b.RAW").write_text("raw")
    write_file(root / "sub1" / "deep" / "c.mzML", seed=3, drift=-4)
    conv = base / "conv.py"   # stand-in converter: -i=<raw> -o=<dir> -> copies the prepared mzML
    conv.write_text("import sys, pathlib, shutil\ni = pathlib.Path(sys.argv[1][3:]); o = pathlib.Path(sys.argv[2][3:])\n"
                    f"shutil.copyfile(pathlib.Path(r'{bank}') / (i.stem + '.mzML'), o / (i.stem + '.mzML'))\n")
    settings = fo.load_settings(root)
    settings["vendor"] = "thermo"
    settings["convert"]["thermo_parser"] = f'"{sys.executable}" "{conv}"'
    fo.save_settings(root, settings)
    for folder, text in ((root, COMPONENTS), (root / "sub1" / "deep", ONLY_CAFFEINE)):
        (folder / fo.SUB).mkdir(exist_ok=True)
        (folder / fo.SUB / "components.yaml").write_text(text)
    return root


def test_units_found(tree):
    units = fo.find_units(tree)
    assert [u.relative_to(tree).as_posix() for u in units] == [".", "sub1", "sub1/deep"]   # not 'empty', not 'MSparser'
    assert fo.find_units(tree, include_subfolders=False) == [tree]
    assert fo.label(tree / "sub1" / "deep", tree) == str(Path("root") / "sub1" / "deep") and fo.label(tree, tree) == "root"
    s = fo.unit_status(tree / "sub1")
    assert (s["raw_files"], s["usable_files"], s["compounds"]) == (1, 0, 0)


def test_run_all_keeps_everything_in_each_folder(tree):
    units = fo.find_units(tree)
    res = fo.run_all(tree, units, fo.load_settings(tree)).set_index("data_folder")
    sub1, deep = str(Path("root") / "sub1"), str(Path("root") / "sub1" / "deep")
    # sub1 has no compound list yet: skipped, the others are done
    assert res.loc["root", "status"] == "done" and res.loc[deep, "status"] == "done"
    assert res.loc[sub1, "status"] == "skipped" and "compound list" in res.loc[sub1, "message"]
    # converted file, settings, compound list, ions and results inside the folder of the raw file
    m = tree / fo.SUB
    for f in ("mzml/a.mzML", "config.yaml", "components.yaml", "work/targets.csv",
              "results/intensity_matrix.csv", "results/intensity_values.csv"):
        assert (m / f).exists(), f
    assert not (tree / "a.mzML").exists()                       # nothing written beside the raw files
    assert (tree / "sub1" / "deep" / fo.SUB / "results" / "intensity_values.csv").exists()
    v = pd.read_csv(m / "results" / "intensity_values.csv").set_index("target_id")
    assert v.loc[CAF, "intensity"] == pytest.approx(1e7, rel=0.03)

    # copy the compound list of root to sub1, then run again: only sub1 is new work
    assert fo.copy_components(tree, [tree / "sub1"]) == 1
    res = fo.run_all(tree, units, fo.load_settings(tree)).set_index("data_folder")
    assert (res["status"] == "done").all()
    assert (tree / "sub1" / fo.SUB / "mzml" / "b.mzML").exists()
    assert "0 processed, 1 unchanged (reused)" in res.loc["root", "log"]
    vb = pd.read_csv(tree / "sub1" / fo.SUB / "results" / "intensity_values.csv").set_index("target_id")
    assert vb.loc[CAF, "mz_window_ppm"] == 10                    # 8 ppm off: found with the second window

    # combined tables at the top folder
    allv = pd.read_csv(m / "results_all" / "intensity_values_all.csv")
    assert sorted(allv["data_folder"].unique()) == sorted(["root", sub1, deep])
    mx = pd.read_csv(m / "results_all" / "intensity_matrix_all.csv").set_index("data_folder")
    assert len(mx) == 3 and mx.loc["root", CAF] == pytest.approx(1e7, rel=0.03)
    # each folder has its own compound list: 'deep' only follows caffeine
    assert pd.isna(mx.loc[deep, "choline [M]+"]) and mx.loc["root", "choline [M]+"] > 0
    assert set(allv.loc[allv["data_folder"] == deep, "component"].dropna()) == {"caffeine"}


def test_ions_follow_the_compound_list(tree):
    deep = tree / "sub1" / "deep"
    settings = fo.load_settings(tree)
    before = pd.read_csv(deep / fo.SUB / "work" / "targets.csv")
    assert set(before["component"]) == {"caffeine"}
    fo.copy_components(tree, [deep])
    fo.run_unit(deep, settings)
    after = pd.read_csv(deep / fo.SUB / "work" / "targets.csv")
    assert "choline" in set(after["component"])
    (deep / fo.SUB / "components.yaml").write_text(ONLY_CAFFEINE)     # back, for the other tests
    fo.run_unit(deep, settings)


def test_a_folder_can_be_processed_alone(tree, capsys):
    """<folder>/MSparser/config.yaml is a normal configuration."""
    cfg = yaml.safe_load((tree / "sub1" / fo.SUB / "config.yaml").read_text())
    assert cfg["paths"]["raw_dir"] == ".." and cfg["paths"]["mzml_dir"] == "mzml"
    main(["intensity", "-c", str(tree / "sub1" / fo.SUB / "config.yaml")])
    assert "1 file(s)" in capsys.readouterr().out


def test_folder_command(tree, capsys):
    main(["folder", str(tree), "--no-convert"])
    out = capsys.readouterr().out
    assert "root: done" in out and "combined tables" in out
    main(["folder", str(tree / "sub1"), "--no-subfolders", "--no-convert"])
    assert "sub1: done" in capsys.readouterr().out


def test_interface(tree):
    AppTest = pytest.importorskip("streamlit.testing.v1").AppTest
    app = Path(__file__).parent.parent / "msparser" / "app.py"
    old = sys.argv
    sys.argv = ["app.py", str(tree), "intensity"]
    try:
        at = AppTest.from_file(str(app), default_timeout=180).run()
        assert not at.exception, [e.value for e in at.exception]
        assert [t.label for t in at.tabs] == ["1. Folders", "2. Compounds", "3. Ions", "4. Convert and extract",
                                              "5. Results"]
        at = [b for b in at.button if b.label == "Extract only"][0].click().run()
        assert not at.exception and not at.error, [e.value for e in at.error]
        assert any("3 of 3 folder(s) done" in s.value for s in at.success)
        labels = {m.label: m.value for m in at.metric}
        assert labels["Folders"] == "3" and labels["Files"] == "3"
        # one folder: its own tables, choice of peaks and detail
        at = [s for s in at.selectbox if s.label == "Results of"][0].set_value("root").run()
        assert not at.exception, [e.value for e in at.exception]
        at = [b for b in at.button if b.label == "Apply the choice of peaks"][0].click().run()
        assert not at.exception and not at.error, [e.value for e in at.error]
    finally:
        sys.argv = old


def test_detected_mz_and_mass_error_are_reported(tree):
    res = tree / "sub1" / fo.SUB / "results"            # the file with masses 8 ppm too high
    v = pd.read_csv(res / "intensity_values.csv").set_index("target_id").loc[CAF]
    assert v["mz_theoretical"] == pytest.approx(195.08765, abs=2e-5)
    assert v["mz_observed"] == pytest.approx(195.08765 * (1 + 8e-6), abs=1e-4)
    assert v["ppm_error"] == pytest.approx(8, abs=0.5)
    cols = list(pd.read_csv(res / "intensity_values.csv").columns)
    assert cols.index("mz_theoretical") < cols.index("mz_observed") < cols.index("ppm_error") < cols.index("flags")
    assert pd.read_csv(res / "intensity_matrix_mz_detected.csv").set_index("file").loc["b", CAF] == pytest.approx(v["mz_observed"])
    assert pd.read_csv(res / "intensity_matrix_ppm_error.csv").set_index("file").loc["b", CAF] == pytest.approx(v["ppm_error"])
    ions = pd.read_csv(res / "intensity_ions.csv").set_index("target_id")
    assert ions.loc[CAF, "ppm_error_median"] == pytest.approx(8, abs=0.5) and ions.loc[CAF, "files_detected"] == 1
    assert ions.loc["caffeine [M+2H]2+", "files_detected"] == 0
    all_ = tree / fo.SUB / "results_all"
    mz = pd.read_csv(all_ / "intensity_matrix_mz_detected_all.csv").set_index("file")
    ppm = pd.read_csv(all_ / "intensity_matrix_ppm_error_all.csv").set_index("file")
    assert mz.loc["a", CAF] == pytest.approx(195.08765, abs=3e-4) and abs(ppm.loc["a", CAF]) < 1


def test_charge_choices():
    from msparser.ui_shared import expand_charges

    assert expand_charges(3, 0, "proton") == [1, 2, 3]          # a molecule takes the lower charges too
    assert expand_charges(1, -2, "proton") == [1, -1, -2]
    assert expand_charges(0, 0, "proton") == []
    assert expand_charges(2, 0, "intrinsic") == [2]             # an ion has exactly its charge
