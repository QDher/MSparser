"""Entry point of the packaged application: starts the interface in the browser.

Usage: mzstab-app [project folder]   (default: <home>/mzstab_project)
"""
import sys
from pathlib import Path


def main():
    base = Path(getattr(sys, "_MEIPASS", Path(__file__).resolve().parent.parent))
    app = base / "mzstab" / "app.py"
    project = sys.argv[1] if len(sys.argv) > 1 else str(Path.home() / "mzstab_project")
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
