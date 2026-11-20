"""External tools; no solver binary or machine-specific path is bundled."""
from pathlib import Path
import os
ROOT=Path(__file__).resolve().parent
BIN=Path(os.environ.get('DUALSPHYSICS_BIN',str(ROOT.parent/'external/DualSPHysics/bin/windows'))).expanduser().resolve()
SOLVER=Path(os.environ.get('DUALSPHYSICS_CHRONO_EXE',str(BIN/'DualSPHysics5.4_chrono.exe'))).expanduser().resolve()
