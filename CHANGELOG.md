# Modification history

## 2026-11-20 — Gunnerus focused-wave case and modelling notes

- Added completed case 05: 16 model seconds, 9,985,001 particles, 0.05 m spacing, free six-DOF vessel with assumed world-fixed 30 N CG force. Added its independent 2-D incident-wave check and measured wave gauges.
- Included explicit reference XML, piston signal, physical inputs and selected sealed-hull/deckhouse contact meshes derived from the MIT vessel.js visual model, with the upstream MIT notice. Turtle geometry remains excluded. Original Gunnerus GLB is external and hash-identified.
- Added a separate `gunnerus/` setup/monitor/render workflow, retaining the earlier turtle scripts unchanged. Adapted external tool/font paths, added a solver-output overwrite guard, extracted unchanged XML helpers and added a reference-input stager. No governing-equation or GPU-kernel edits.
- Published both accepted 26 s portrait and landscape films, actual motion/load CSV, run audits, resolution and attitude summaries. Excluded full particle fields, render caches, solver binaries, environments and exploratory cases.
- Added notes on rigid-body mass/CG/inertia, Froude versus Reynolds scaling, SPH/FVM/IBM/DEM, time/particle resolution and possible VRes or RANS-map extensions. VRes, actuator disk, RANS coupling, flooding and capsize validation are not implemented in this result.
- Packaging verification is documented separately from the completed simulation; no new CFD study was launched for this publication.

## 2026-11-18 — public demonstration package

- Added source-specific local hull reconstruction, history-informed mass/CG/inertia assumptions and hydrostatic equilibrium preparation.
- Added Capytaine 3 six-DOF radiation/diffraction calculation, response assembly and evaluated surface visualization.
- Added DualSPHysics piston-tank inputs, DBC walls/piston, mDBC hull normals and initially moving/free six-DOF hull. Native solver equations and kernels are unchanged.
- Added wider focused-wave case after empty-flume commissioning; corrected floor/piston particle ownership and extended the stationary floor behind the piston stroke.
- Added particle-based rigid motion extraction, density/velocity/roof-contact diagnostics, smoothed water surfaces and short English portrait story.
- Published the accepted 35-second film with its measured three-second bow-impact intro, input snapshots and compact results. Excluded failed prototypes, full particle fields, solver binaries, private source/derived geometry and credentials.
- Adapted resource paths for `TURTLE_ASSET_DIR`, `DUALSPHYSICS_BIN` and `SHIP_FONT_DIR`; sanitized reference metadata and XML placeholders. Added handoff documentation. Full CFD and rendering were not rerun just for publication.
