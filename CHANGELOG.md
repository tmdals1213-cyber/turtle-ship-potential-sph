# Modification history

## 2026-11-18 — public demonstration package

- Added source-specific local hull reconstruction, history-informed mass/CG/inertia assumptions and hydrostatic equilibrium preparation.
- Added Capytaine 3 six-DOF radiation/diffraction calculation, response assembly and evaluated surface visualization.
- Added DualSPHysics piston-tank inputs, DBC walls/piston, mDBC hull normals and initially moving/free six-DOF hull. Native solver equations and kernels are unchanged.
- Added wider focused-wave case after empty-flume commissioning; corrected floor/piston particle ownership and extended the stationary floor behind the piston stroke.
- Added particle-based rigid motion extraction, density/velocity/roof-contact diagnostics, smoothed water surfaces and short English portrait story.
- Published the accepted 35-second film with its measured three-second bow-impact intro, input snapshots and compact results. Excluded failed prototypes, full particle fields, solver binaries, private source/derived geometry and credentials.
- Adapted resource paths for `TURTLE_ASSET_DIR`, `DUALSPHYSICS_BIN` and `SHIP_FONT_DIR`; sanitized reference metadata and XML placeholders. Added handoff documentation. Full CFD and rendering were not rerun just for publication.
