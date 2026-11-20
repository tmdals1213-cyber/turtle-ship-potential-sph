# Gunnerus workflow scripts

This directory contains the Python setup, monitoring and rendering workflow used for the completed focused-wave example. The numerical solver is **upstream DualSPHysics**, not a new solver implemented here. The earlier turtle workflow in `study/` is preserved.

Start with [the case](../SPH_CFD-incomp-free-surface-6DOF_GPU_Cpp-CUDA-Python_05-gunnerus-focused-wave), [reproduction](../docs/GUNNERUS_REPRODUCE.md) and [physics notes](../docs/MARINE_MODEL_NOTES.md). `prepare_reference.py` stages the accepted inputs without running a simulation. Set `PYTHONUTF8=1` on Windows.

| Files | Purpose |
|---|---|
| `prepare_reference.py` | Install reference geometry and inputs into ignored local working directories |
| `prepare_gunnerus*.py`, `place_gunnerus_contact.py`, `xml_helpers.py` | Geometry/loading and tank/wavemaker generation; optional variants |
| `run_sph.py`, `check_gunnerus_input.py`, `check_gunnerus_focus_wave.py` | External solver launch and input/wave checks |
| `read_native_float.py`, `sph_segments.py`, `analyze_propulsion.py`, `report_gunnerus*.py` | Native rigid-body records, forces, pose and monitors |
| `prepare_gunnerus_surface_pressure.py` | Invoke upstream MeasureTool for nearby-water pressure |
| `render_gunnerus*.py`, `add_method_comparison.py`, `finish_gunnerus_focus.py` | Actual-data replay, portrait/landscape composition and decoding checks |

Several modules retain functions for the earlier calm-water pilots. Use the final-case commands in the reproduction document; running every script with defaults would select those older pilots. Native restart arguments are retained for provenance, but the published final case ran uninterrupted. Restart repair tools and unsuccessful prototypes are not part of the supported reproduction path.

[SOURCE_ADAPTATIONS.json](SOURCE_ADAPTATIONS.json) records source hashes and publication changes: external tool/font paths, an overwrite guard and extraction of two unchanged XML helpers. Physics inputs and numerical equations were not changed for publication. The reference staging helper and documentation are new packaging files.
