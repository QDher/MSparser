# PyInstaller specification. Build from the repository root:
#   pyinstaller --noconfirm packaging/mzstab.spec
# The result is the folder dist/mzstab-app (executable + libraries).
# An executable only runs on the operating system it was built on.
from pathlib import Path

from PyInstaller.utils.hooks import collect_all, collect_submodules, copy_metadata

root = Path(SPECPATH).parent
datas = [(str(root / "mzstab" / "app.py"), "mzstab")]
binaries, hidden = [], collect_submodules("mzstab")
for pkg in ("streamlit", "altair", "pyarrow", "ruptures", "molmass", "pyteomics", "psims", "statsmodels", "patsy"):
    d, b, h = collect_all(pkg)
    datas += d
    binaries += b
    hidden += h
for pkg in ("streamlit", "altair", "pandas", "numpy", "pyarrow"):
    datas += copy_metadata(pkg)

a = Analysis(
    [str(root / "packaging" / "run_app.py")],
    pathex=[str(root)], binaries=binaries, datas=datas, hiddenimports=hidden,
    excludes=["pytest", "IPython", "notebook", "tkinter"],
)
pyz = PYZ(a.pure)
exe = EXE(pyz, a.scripts, [], exclude_binaries=True, name="mzstab-app", console=True)
coll = COLLECT(exe, a.binaries, a.datas, name="mzstab-app")
