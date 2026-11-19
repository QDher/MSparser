"""Entry point of the packaged application: starts the interface in the browser.

Usage: MSparser [project folder]   (default: <home>/MSparser_project)
"""
import sys
from pathlib import Path


def main():
    base = Path(getattr(sys, "_MEIPASS", Path(__file__).resolve().parent.parent))
    app = base / "msparser" / "app.py"
    default, former = Path.home() / "MSparser_project", Path.home() / "mzstab_project"
    if not default.exists() and former.exists():
        default = former  # project folder created by a version named mzstab
    project = sys.argv[1] if len(sys.argv) > 1 else str(default)
    from streamlit.web import cli as stcli

    sys.argv = [
        "streamlit", "run", str(app),
        "--global.developmentMode=false",
        "--browser.gatherUsageStats=false",
        "--server.fileWatcherType=none",
        "--", project,
    ]
    sys.exit(stcli.main())


if __name__ == "__main__":
    main()
