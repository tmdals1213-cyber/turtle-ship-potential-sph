# Reproducing the examples

Run commands from the repository root. Reference environment: Windows, Python 3.13, Ryzen 9 9950X, RTX 5080, CUDA 12.8. The scripts and inputs are supplied; raw fields and third-party model files are not. Running the calculations is separate from downloading the finished videos.

## 1. Python and external solver

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
python -m pip install -r requirements-render.txt
```

Install FFmpeg separately and put it on PATH for video encoding. The reference renderers use Windows Segoe UI, Segoe UI Bold, Segoe UI Symbol and Malgun fonts. `SHIP_FONT_DIR` can point to a directory containing those font filenames. Non-Windows solver launch and font substitution are not verified.

Obtain [DualSPHysics](https://github.com/DualSPHysics/DualSPHysics) separately. Pinned source revision: `ef3721a861fda961f0e2f9ec4cd317b19de99086`. The reference solver reported 5.4.355 and GenCase 5.4.354.01. Set:

```powershell
$env:DUALSPHYSICS_BIN = 'D:\tools\DualSPHysics\bin\windows'
```

The launcher expects `GenCase_win64.exe` and `DualSPHysics5.4_win64.exe`, with their runtime dependencies in that directory. The recorded RTX 5080 run used a local CUDA 12.8 / SM 120 / MSVC 2022 build, C++14 and CUDA14, static CUDA runtime, OpenMP, and `--use_fast_math`. Build definitions included `_WITHMR`, `_WITHGPU`, `DISABLE_CHRONO`, `DISABLE_MOORDYNPLUS`, `DISABLE_WAVEGEN`, `_CRT_SECURE_NO_WARNINGS`, `CUB_IGNORE_DEPRECATED_CPP_DIALECT` and `THRUST_IGNORE_DEPRECATED_CPP_DIALECT`. Only the build configuration was adapted; governing equations and numerical kernels were unchanged. The optional internal WaveGen library is not required because piston motion is supplied as a time series. The upstream binaries may require a different GPU architecture; a tested SM 120 binary is not included here.

## 2. Asset-free commissioning

```powershell
python study/prepare_focused_turtle.py --name focus_empty_A065_near_40 --empty-2d --amplitude 0.65 --dp 0.04 --duration 12 --focus-x 1 --focus-time 7.5 --peak-period 2.2 --piston-x 6
python study/run_sph.py focus_empty_A065_near_40 --generate-only
python study/run_sph.py focus_empty_A065_near_40 --solver-only --guard
python study/analyze_focus_empty.py focus_empty_A065_near_40
```

Check the generated geometry and gauges before adding a hull. The published empty-flume result is a reference, not a substitute for checking a changed wavemaker. The input target `0.65 m` is a **linear coherent crest amplitude**, not measured wave height or significant wave height.

## 3. Prepare the private turtle geometry

Obtain the Blender archive from the [creator](https://wifewell.tistory.com/entry/거북선-3D파일FBXOBJ-및-블렌더-무료-나눔); source/derived geometry must remain local. Default asset directory: `private_assets/turtle/`, ignored by Git. Alternatively set `TURTLE_ASSET_DIR` to an absolute directory outside this repository.

The original asset uses Blender 2.79. Convert it using a compatible Blender, with embedded-script auto-execution disabled. Replace the example paths below with your own paths:

```powershell
$env:TURTLE_ASSET_DIR = Join-Path (Get-Location) 'private_assets\turtle'
& 'D:\tools\blender279\blender.exe' --disable-autoexec --background 'D:\assets\turtle.blend' --python study/export_turtle_blender.py -- $env:TURTLE_ASSET_DIR --export
python study/prepare_turtle_visual.py
python study/reconstruct_turtle.py --extended
python study/prepare_history_model.py
python study/prepare_turtle_materials.py
```

These scripts are specific to the object names and coordinates of that source asset. They do not accept arbitrary ship CAD as a drop-in replacement. If adapting a different hull, recompute buoyancy, mass placement, inertia, normals and BEM panels. Keep all sampled half-breadths, panels, closed STL, source meshes, roof maps and render caches private too. Saving the downloaded archive as `turtle_blend.7z` in the asset directory allows its optional SHA-256 provenance check.

## 4. Potential-flow result

```powershell
python study/run_history_potential.py
```

Outputs are under `study/results/potential_turtle_history/`. The BEM uses the mean wetted hull, six harmonic DOFs and fixed reference speed. The plotted incident surface includes a second-order Stokes shape plus linear BEM scattering/radiation; this is an approximate visualization, not a fully nonlinear potential solution. The public `01-potential` folder preserves the original summary separately.

## 5. Coasting SPH

```powershell
python study/prepare_tow.py --name free6_history_40 --history --free6 --mixed-walls --initial-speed 0.5 --height 0.24 --dp 0.04 --duration 16 --coefh 1.2 --slip-mode 1
python study/run_sph.py free6_history_40 --generate-only
python study/audit_history_particles.py free6_history_40
python study/run_sph.py free6_history_40 --solver-only --guard
python study/summarize_sph.py free6_history_40
python study/process_free6.py free6_history_40
python study/report_free6.py
```

Do not add `--fitted`: `--history` selects the final physical model itself. The initial velocity is +x, waves propagate toward -x; there is no continuing thrust or prescribed surge. The native rigid-body solver handles all six DOFs. Preserve prior output when changing inputs. The launcher refuses to replace an existing solver-output directory.

## 6. Focused-wave SPH

```powershell
python study/prepare_focused_turtle.py --name focus_turtle_A065_near_40 --amplitude 0.65 --dp 0.04 --duration 10 --focus-x 1 --focus-time 7.5 --peak-period 2.2 --piston-x 6
python study/run_sph.py focus_turtle_A065_near_40 --generate-only
python study/audit_history_particles.py focus_turtle_A065_near_40
python study/run_sph.py focus_turtle_A065_near_40 --solver-only --guard
python study/process_focused_turtle.py focus_turtle_A065_near_40
python study/report_focused_turtle.py focus_turtle_A065_near_40
```

The stationary floor and side walls extend behind the piston stroke, with static-wall particle ownership at shared points. This avoids the leakage found in an earlier rejected commissioning setup. Do not silently retune the physical mass to remove the coarse SPH equilibrium bias.

## 7. Re-render after producing local fields

```powershell
python study/prepare_portrait_scene.py
python study/render_turtle_v6_portrait.py
python study/render_focused_turtle.py focus_turtle_A065_near_40 --start 4.60 --end 5.92 --seconds 3 --name turtle_focused_intro
python study/assemble_focused_story.py
```

These require the private visual model, regenerated BEM surface/panels, both SPH particle sequences and local render caches. They cannot reconstruct fluid fields from the compact published motion JSON alone. Output is under ignored `study/media/`; accepted public films remain in root `media/`. Full simulation and full rendering have not been rerun solely for packaging.

## Reference input snapshots

The top-level case folders contain the recorded parameters and piston signals. Their XML files use `@HULL_STL@` and `@HULL_NORMALS@` placeholders deliberately: restricted geometry is not embedded. Use the generators above to create local absolute hull paths and correctly relative normal paths. Compare generated physical settings with the reference snapshots before an expensive run; do not feed unresolved placeholders to GenCase.
