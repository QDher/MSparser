#!/usr/bin/env bash
# Build the Linux or macOS application. Run from the repository root:  bash packaging/build.sh
# Result: dist/mzstab-app/mzstab-app
set -euo pipefail
python3 -m venv .build-venv
source .build-venv/bin/activate
python -m pip install --upgrade pip
pip install ".[ui]" pyinstaller
pyinstaller --noconfirm packaging/mzstab.spec
echo "Built: dist/mzstab-app/mzstab-app  (distribute the whole dist/mzstab-app folder)"
