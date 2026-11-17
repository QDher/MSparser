"""Synthetic centroided mzML files for testing (requires psims)."""
from __future__ import annotations

from pathlib import Path

import numpy as np
from psims.mzml.writer import MzMLWriter

# (m/z, plateau intensity, isotope shift, isotope ratio)
C13 = 1.00335483534
IONS = {
    "caffeine_H": (195.087652, 1e7, C13, 8 * 0.0107 / 0.9893),
    "caffeine_Na": (217.069597, 1e6, C13, 8 * 0.0107 / 0.9893),
    "choline": (104.106990, 5e6, C13, 5 * 0.0107 / 0.9893),
    "glucose_Na": (203.052608, 2e5, C13, 6 * 0.0107 / 0.9893),
}


def profile(t, drift=0.0):
    """Rise 0-5 s, plateau 5-25 s, decay afterwards; optional linear drift (%/s)."""
    p = np.clip(t / 5.0, 0, 1) * np.clip((30.0 - t) / 5.0, 0, 1)
    return p * (1 + drift / 100 * (t - 15))


def write_file(path: Path, seed: int, drift: float = 0.0, interference: bool = False,
               ions=IONS, n_scans: int = 90, duration: float = 30.0):
    rng = np.random.default_rng(seed)
    times = np.linspace(0.2, duration - 0.2, n_scans)
    with MzMLWriter(open(path, "wb"), close=True) as out:
        out.controlled_vocabularies()
        out.file_description(["MS1 spectrum"])
        out.software_list([{"id": "psims-writer", "version": "0.1", "params": ["python-psims"]}])
        src = out.Source(1, ["electrospray ionization"])
        ana = out.Analyzer(2, ["orbitrap"])
        det = out.Detector(3, ["inductive detector"])
        ic = out.InstrumentConfiguration(id="IC1", component_list=[src, ana, det])
        out.instrument_configuration_list([ic])
        pm = out.ProcessingMethod(order=1, software_reference="psims-writer",
                                  params=["Conversion to mzML"])
        out.data_processing_list([out.DataProcessing([pm], id="DP1")])
        with out.run(id="run1", instrument_configuration="IC1"):
            with out.spectrum_list(count=n_scans + 2):
                for i, t in enumerate(times):
                    mz, inten = [], []
                    for m, top, shift, ratio in ions.values():
                        level = top * profile(t, drift)
                        for mm, lv in ((m, level), (m + shift, level * ratio)):
                            if lv < 2e3:  # noise threshold: peak not reported
                                continue
                            rsd = np.sqrt(2e3 / lv + 0.02**2)
                            mz.append(mm * (1 + rng.normal(0, 1e-6)))
                            inten.append(max(lv * (1 + rng.normal(0, rsd)), 1.0))
                    if interference:  # isobaric peak 7 ppm above caffeine [M+H]+
                        mz.append(195.087652 * (1 + 7e-6))
                        inten.append(3e5 * (1 + rng.normal(0, 0.05)))
                    for _ in range(30):  # background
                        mz.append(rng.uniform(60, 740))
                        inten.append(rng.uniform(2e3, 6e3))
                    order = np.argsort(mz)
                    mz, inten = np.array(mz)[order], np.array(inten)[order]
                    out.write_spectrum(
                        mz, inten, id=f"controllerType=0 controllerNumber=1 scan={i + 1}",
                        params=[{"ms level": 1}, {"total ion current": float(inten.sum())}],
                        centroided=True, polarity=1, scan_start_time=t / 60.0,
                        scan_window_list=[(50.0, 750.0)],
                        scan_params=[{"ion injection time": float(rng.uniform(5, 20))},
                                     {"mass resolving power": 140000.0}],
                    )
                for j in range(2):  # MS2 scans, to be dropped by the scan filter
                    out.write_spectrum(
                        np.array([100.0, 150.0]), np.array([1e4, 2e4]),
                        id=f"controllerType=0 controllerNumber=1 scan={n_scans + j + 1}",
                        params=[{"ms level": 2}], centroided=True, polarity=1,
                        scan_start_time=(duration + j) / 60.0, scan_window_list=[(50.0, 750.0)],
                    )


COMPONENTS = """\
components:
  - {name: caffeine, formula: C8H10N4O2, charges: [1, 2], mode: proton}
  - {name: glucose, formula: C6H12O6, charges: [1], mode: proton}
  - {name: sodium, formula: Na, charges: [1], mode: intrinsic}
  - {name: choline, formula: C5H14NO, charges: [1], mode: intrinsic}
"""

CONFIG = """\
vendor: thermo
extraction:
  ppm_windows: [10, 5]
resolution: {value: 60000, at_mz: 200, scaling: orbitrap}
thresholds: {rsd_factor: auto, slope_limit: auto}
windows:
  lengths_s: [3, 5, 8, 10, 15]
  step_s: 1.0
  min_scans: 8
"""


def make_project(root: Path, n_files: int = 8) -> Path:
    (root / "mzml").mkdir(parents=True, exist_ok=True)
    rows = ["file,reaction,components,run_order"]
    for i in range(n_files):
        drift = -4.0 if i == n_files - 1 else 0.0  # last file drifts throughout
        write_file(root / "mzml" / f"run{i + 1:02d}.mzML", seed=i, drift=drift,
                   interference=(i % 4 == 0))
        rows.append(f"run{i + 1:02d},{'A' if i % 2 else 'B'},,{i + 1}")
    (root / "samples.csv").write_text("\n".join(rows) + "\n")
    (root / "components.yaml").write_text(COMPONENTS)
    (root / "config.yaml").write_text(CONFIG)
    return root / "config.yaml"


if __name__ == "__main__":
    import sys

    print(make_project(Path(sys.argv[1])))
