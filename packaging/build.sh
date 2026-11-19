#!/usr/bin/env bash
# Build the Linux or macOS application. Run from the repository root:  bash packaging/build.sh
# Result: dist/MSparser/MSparser
set -euo pipefail
python3 -m venv .build-venv
source .build-venv/bin/activate
python -m pip install --upgrade pip
pip install ".[ui]" pyinstaller
pyinstaller --noconfirm packaging/msparser.spec
echo "Built: dist/MSparser/MSparser  (distribute the whole dist/MSparser folder)"
