# Reproducing the tanker example

Run from the repository root. The reference platform was Windows, Python 3.13, CUDA 12.8, RTX 5080 and MSVC 2022. No third-party geometry is needed. Python dependencies are in the root requirements files; FFmpeg must be on PATH. `TANKER_FONT_DIR` may point to a directory containing `arial.ttf` and `arialbd.ttf`. Non-Windows launch/font substitutions are not verified.

## External software

Obtain the pinned DualSPHysics revision in THIRD_PARTY_NOTICES. Set a directory containing `GenCase_win64.exe`, `PartVTK_win64.exe`, and `ComputeForces_win64.exe`. Obtain ComputeForces from the [official tool directory](https://dual.sphysics.org/sphcourse/DualSPHysics-bin/dualsphysics/bin/). Its recorded version is 5.4.266.02, SHA-256 `59d6a56c6ea1cb591f7015187b9fd5c06b47122192bb33e5c6f8ebbb940e5c2c`. It is a separate unmodified freeware postprocessor, not a source component relicensed by this repository.

The tanker requires a **Chrono-enabled GPU solver**, unlike the older ship executable. Set its full executable path explicitly:

```powershell
$env:DUALSPHYSICS_BIN = 'D:\tools\DualSPHysics\bin\windows'
$env:DUALSPHYSICS_CHRONO_EXE = 'D:\tools\tanker-chrono\DualSPHysics5.4_chrono.exe'
python -m pip install -r requirements.txt -r requirements-render.txt
```

The solver directory must also contain matching `dsphchrono.dll`, `ChronoEngine.dll`, OpenMP runtime and `DsphConfig.xml`. They remain external. The reference build enabled native Chrono and disabled only MoorDynPlus/internal WaveGen; numerical source was unchanged.

For an equivalent Windows build, use `tanker/build/CMakeLists.txt` from an x64 MSVC developer shell with CUDA 12.8 and CMake/Ninja available:

```powershell
cmake -S tanker/build -B external/tanker-build -G Ninja -DDSPH_SOURCE=D:/tools/DualSPHysics/src/source -DCMAKE_BUILD_TYPE=Release -DCMAKE_CUDA_ARCHITECTURES=120
cmake --build external/tanker-build --parallel 6
```

This recipe expects the upstream `src/lib/vs2022/dsphchrono.lib` matching the runtime DLLs. Adjust CUDA architecture to the target GPU. It preserves the recorded source list and compiler definitions; toolchain installation is not automated. The reference build took 483 s, but the portable recipe was not rebuilt during packaging.

## Generate, inspect, then run

```powershell
python -X utf8 tanker/prepare_reference.py
```

This creates procedural geometry, source particle layouts and all four v8 vehicle inputs under ignored `tanker/geometry/` and `tanker/sph/`. It runs GenCase and particle audits, **not CFD**. Each fluid count must be 137,812. Fixed source-layout folders are geometry stages; their fixed-tank pilot conditions are not the accepted coupled comparison.

```powershell
python -X utf8 tanker/run_reference.py
```

Runs all four accepted inputs sequentially, checks completion/particle exclusion/rigid pose/partition integrity, and extracts grouped loads. Existing output is preserved. A failed case is reported while the remaining cases continue. Existing run records do not imply raw fields are still present; after cleanup, use a clean working clone/folder before regeneration. Do not remove immutable top-level reference snapshots.

The launch guard stops exclusions exceeding 0.1% or persistent minimum-time-step clamps. These are failure guards, not a complete validation criterion. Inspect `fluid_boundary_qa.json`, `coupled_analysis.json` and `surface_loads/analysis.json`. The summary uses native forces and the same dry-vehicle inertia/CG across cases. Failed pilots and interrupted runs are not accepted comparison data.

## Re-render

```powershell
python -X utf8 tanker/render_final_edition.py portrait
python -X utf8 tanker/render_final_edition.py landscape
```

Both need regenerated full particle fields and loads. They produce ignored `tanker/media/` outputs; public final films are in root `media/`. Names and cameras reproduce the accepted v8 edit. The title/conclusion asserts the recorded reference outcome; a changed physical case needs its own name and a revised narrative, never forced motion or favorable results.

## Packaging verification

The publishing pass compiles/imports the portable workflow and regenerates/audits all four input sets against compact reference parameters. It does not repeat the four CFD runs or claim cross-platform/time-resolution validation. Refer to `case 06/package_verification.json` for the actual checks performed.
