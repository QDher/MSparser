"""Entry point of the packaged application: starts the interface in the browser.

Usage: MSparser [folder [intensity|stability]]
"""
import sys
from pathlib import Path


def main():
    base = Path(getattr(sys, "_MEIPASS", Path(__file__).resolve().parent.parent))
    app = base / "msparser" / "app.py"
    extra = sys.argv[1:3]  # optional: a folder to open, and the mode
    from streamlit.web import cli as stcli

    sys.argv = [
        "streamlit", "run", str(app),
        "--global.developmentMode=false",
        "--browser.gatherUsageStats=false",
        "--server.fileWatcherType=none",
        "--", *extra,
    ]
    sys.exit(stcli.main())


if __name__ == "__main__":
    main()
