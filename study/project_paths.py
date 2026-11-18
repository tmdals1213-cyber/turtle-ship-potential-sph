# SPDX-License-Identifier: GPL-3.0-or-later
"""Local resource locations; geometry and generated fields are never published."""
from pathlib import Path
import os

STUDY = Path(__file__).resolve().parent
REPO = STUDY.parent
ASSET_DIR = Path(os.environ.get('TURTLE_ASSET_DIR', REPO / 'private_assets/turtle')).expanduser().resolve()
SOLVER_BIN = Path(os.environ.get('DUALSPHYSICS_BIN', REPO / 'external/DualSPHysics/bin/windows')).expanduser().resolve()
FONT_DIR = Path(os.environ.get('SHIP_FONT_DIR', 'C:/Windows/Fonts')).expanduser()
for directory in [ASSET_DIR, STUDY / 'results', STUDY / 'sph', STUDY / 'media']:
    directory.mkdir(parents=True, exist_ok=True)
