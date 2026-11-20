# Reproduce the Gunnerus focused-wave example

Run from the repository root on Windows. The accepted [case snapshots](../SPH_CFD-incomp-free-surface-6DOF_GPU_Cpp-CUDA-Python_05-gunnerus-focused-wave) and videos are preserved separately from newly generated outputs. Do not launch an expensive run just to read the results.

## Dependencies

The same upstream DualSPHysics revision/build described in [REPRODUCE.md](REPRODUCE.md#1-python-and-external-solver) was used: source `ef3721a861fda961f0e2f9ec4cd317b19de99086`, solver 5.4.355, CUDA 12.8 / SM 120 on RTX 5080. Governing equations and numerical kernels were not edited. Obtain GenCase, the GPU solver and MeasureTool with their runtime libraries separately. No binaries are bundled.

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt -r requirements-render.txt
python -m pip install manifold3d
$env:PYTHONUTF8 = '1'
$env:DUALSPHYSICS_BIN = 'D:\tools\DualSPHysics\bin\windows'
```

`manifold3d` is only needed to regenerate the optional deckhouse union. Reference staging uses only the Python standard library. FFmpeg must be on PATH for videos. The font override `SHIP_FONT_DIR` should contain `segoeui.ttf` and `segoeuib.ttf`; the default is Windows Fonts. Cross-platform solver launch/rendering is not verified.

## Stage and calculate accepted inputs

```powershell
python gunnerus/prepare_reference.py
$ship = 'gunnerus_focus_A115_x5.5_D300_050_16s'
$empty = 'gunnerus_focus_empty_A115_D300_050_16s'

python gunnerus/run_sph.py $empty --generate-only
python gunnerus/check_gunnerus_input.py $empty
python gunnerus/run_sph.py $empty --solver-only --guard
python gunnerus/check_gunnerus_focus_wave.py $empty

python gunnerus/run_sph.py $ship --generate-only
python gunnerus/check_gunnerus_input.py $ship
python gunnerus/run_sph.py $ship --solver-only --guard
```

`prepare_reference.py` expands `@GUNNERUS_ROOT@` paths and copies the exact MIT-derived reference contact geometry into ignored `gunnerus/geometry/`. It copies **inputs only**, not successful-run status or old results. It refuses existing geometry/case folders. Generated inputs go under ignored `gunnerus/sph/`, leaving the reference snapshots unchanged. Raw particle files can consume substantial disk space and are not included here.

The empty case is a separate 2-D incident-wave check, not the same 3-D ship simulation without a ship. The recorded ship run completed 16.00005 model seconds in about 2 h 52 min; hardware, build and output costs affect rerun time. Inspect the initial audit, `Run.out`, particle losses, time-step clamping and density before accepting a rerun. `--guard` is a failure guard, not an accuracy certificate.

To stop normally, write `0` to the active output directory's `TERMINATE` file and wait for native completion. Keep its output. The accepted run was uninterrupted; restart repair is not part of this reproduction recipe.

## Optional input generation / variants

The generator scripts are included so the inputs are explainable and editable. Reference contact geometry is already supplied; you do not need to regenerate a visual hull to run it. For a distinctly named new wave condition, use:

```powershell
python gunnerus/prepare_gunnerus_focus.py --amplitude 1.10 --station-x 5.5 --period 2.6 --depth 3
```

That is a **new uncomputed case**, not the published result. Existing case directories are rejected. The final reference corresponds to amplitude 1.15, period 2.6 and depth 3; its `setup.json` also records the separately measured wave encounter selection. The target amplitude is a linear coherent crest, not significant wave height.

`prepare_gunnerus.py`, `prepare_gunnerus_impact_geometry.py` and `place_gunnerus_contact.py` document mass/CG/inertia, sealed hull, simplified deckhouses and station placement. The supplied `hull_closed_candidate.stl` is the accepted upstream-visual-derived closed seed. Rebuilding that seed from arbitrary new visual CAD is not automated or claimed; replacing it requires fresh geometry checks. Do not run geometry constructors over existing accepted geometry.

## Visual model for rendering only

The original 31.8 MB GLB is not duplicated in this repository. It is the **2020 vessel.js visual model**, not the newer C-JOB engineering geometry. Download it from the original source and verify the recorded hash before use:

```powershell
$model = 'gunnerus/geometry/candidates/gunnerus_vesseljs/Gunnerus.glb'
Invoke-WebRequest 'https://raw.githubusercontent.com/shiplab/vesseljs/dev/examples/3D_models/GLTF/Gunnerus.glb' -OutFile $model
$expected = '50af42c962c27be109b5f8f24171c801c525e63173feb4d53d1ffbae9a57cf3c'
if ((Get-FileHash $model -Algorithm SHA256).Hash.ToLowerInvariant() -ne $expected) {
    throw 'Upstream model differs from the recorded reference; inspect before rendering.'
}
```

Source Git blob SHA: `0557a1e975646586b588ce0c3374ceee72a7dae8`. See [source notices](../THIRD_PARTY_NOTICES.md#gunnerus-visual-model-and-derived-contact-geometry) and the included MIT notice. Downloading the model is unnecessary for the GPU calculation when using the supplied contact geometry.

## Analyze and render after a successful solver run

```powershell
python gunnerus/finish_gunnerus_focus.py --case $ship --prefix gunnerus_stronger --video-name gunnerus_stronger_waves_story_portrait --dilute-water --surface-level 0.38
python gunnerus/add_method_comparison.py
python gunnerus/render_gunnerus_landscape.py
```

The first command reads native body histories, reconstructs actual particle poses, invokes MeasureTool, builds caches, renders the 22 s precursor and decodes it. The second inserts the four-second conceptual comparison, producing the 26 s portrait. The last re-renders the same result in landscape; it uses intro images and metadata from the previous steps. Inspect the encoded outputs before marking a rerun's visual review passed.

The finished public videos are in root `media/`. New videos, monitors, caches and fields remain ignored under `gunnerus/`. The compact public CSV cannot reconstruct the missing full fluid fields. The exact recorded video bytes are preserved, but new builds/fonts/dependencies can change a rerender.

## Publication verification boundary

Packaging checks cover Python syntax/import availability, staged geometry/XML paths, unchanged physical inputs/piston data and media hashes. The input stager was exercised in an isolated temporary workspace. The full 10-million-particle solver and complete film rendering were **not rerun just for publication**. [SOURCE_ADAPTATIONS.json](../gunnerus/SOURCE_ADAPTATIONS.json) lists the limited portability changes.
