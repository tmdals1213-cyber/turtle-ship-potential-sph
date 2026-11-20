# Sources, credits and redistribution boundaries

## Repository code

Python orchestration, case configurations and documentation are provided under GPL-3.0-or-later; the full text is in LICENSE. No third-party solver source or binary is bundled. The license of this repository does not grant rights to redistribute the external turtle mesh or relicense third-party artwork in the videos.

## DualSPHysics

- Source: <https://github.com/DualSPHysics/DualSPHysics>
- Pinned commit: `ef3721a861fda961f0e2f9ec4cd317b19de99086`
- Downloaded source archive SHA-256: `cd8461cea239e31695a71be94cdd1a573070f4f048e41e814b2ff5f611590362`
- Formulation: <https://github.com/DualSPHysics/DualSPHysics/wiki/3.-SPH-formulation>
- The pinned root LICENSE is LGPL-2.1; source headers include LGPL-2.1-or-later and GPL-3.0-or-later notices. Retain the applicable per-file upstream notices when obtaining/building/distributing that software. This package does not replace them with a single blanket license.
- Reference build: CUDA 12.8, SM 120, MSVC 2022, optional Chrono/MoorDynPlus/internal WaveGen disabled. Only build configuration was adapted. No governing-equation or GPU-kernel edits.

## Capytaine

- Version 3.0.0; Apache-2.0, installed as an external dependency.
- Project: <https://capytaine.org/>
- Theory: <https://capytaine.org/stable/theory_manual/theory.html>
- Used for six-DOF linear radiation/diffraction and harmonic motion at a reference speed. Project-specific changes are the Python geometry, loading, API calls, response assembly and visualization.

## Turtle visual model

Creator: **예은이네 / 와이프잘모시기**. Source: [거북선 3D(FBX,OBJ 및 블렌더 파일) 무료 나눔](https://wifewell.tistory.com/entry/거북선-3D파일FBXOBJ-및-블렌더-무료-나눔), posted 2019-09-15, terms checked 2026-11-18.

The creator permits educational and commercial rendered use, including blogs/YouTube, and prohibits redistribution/sale of the model files. Source and all derived geometry are omitted from both GitHub and transfer archives. Acquire your own source copy from the creator. Credit is retained with the videos and here; no claim of authorship of the third-party model is made.

Source archive SHA-256: `a091c0775f920565258a843a0ca1232206b49234cffa4fca60200a20bb5506c9`. Blender 2.79b was used locally with automatic script execution disabled. [Official Windows archive](https://download.blender.org/release/Blender2.79/blender-2.79b-windows64.zip), SHA-256 `b88d78c7f10bdb532ef511cdf770081eecc983ddc95a08ae96a1f7a77658b55e`.

The public images/videos are rendered outputs, not a downloadable mesh license. Plausible wood/metal colors are artistic choices. The source visual model, fitted symmetric envelope, assumed loading and historical reconstruction reports are distinct sources of information.

## Historical and wave references

- [2026 turtle-ship reconstruction announcement](https://www.kculture.or.kr/brd/board/399/L/menu/498?bbIdx=15970&brdType=R): reported dimensions and self-weight used as scenario inputs, based on 1795 documents rather than proof of the exact 1592 vessel.
- [National Maritime Museum report catalogue](https://www.seamuse.go.kr/resources/academicreport/info/4452): loading analogue and research context; CG/inertia remain declared assumptions.
- [Green water loads on prismatic obstacles, Journal of Fluid Mechanics](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/green-water-loads-on-prismatic-obstacles/DF5D3EB421E38AEF11D70DFD967BEC12): context for focused-wave/green-water studies. That paper's validation is not validation of this turtle hull.

Other libraries are installed separately via the requirements files and retain their respective licenses. Fonts, FFmpeg, Blender and CUDA are not bundled.

## Gunnerus visual model and derived contact geometry

- Source visual model: [vessel.js / Gunnerus.glb](https://github.com/shiplab/vesseljs/blob/dev/examples/3D_models/GLTF/Gunnerus.glb), visual model by IHB NTNU, NTNU Shiplab; added by Felipe Ferrari in 2020.
- Original source license: [MIT, Copyright (c) 2017 shiplab](https://github.com/shiplab/vesseljs/blob/dev/LICENSE.md). The exact notice is included with case 05 geometry as `geometry/candidates/gunnerus_vesseljs/LICENSE.upstream.md` and applies to the derived reference meshes too.
- Original GLB SHA-256: `50af42c962c27be109b5f8f24171c801c525e63173feb4d53d1ffbae9a57cf3c`; Git blob `0557a1e975646586b588ce0c3374ceee72a7dae8`. Original GLB is not bundled.
- [Upstream provenance](https://github.com/shiplab/gunnerus/blob/main/derived/README.md) distinguishes this older approximate vessel.js visual/teaching model from the newer C-JOB engineering geometry with separate terms. The latter was not incorporated here.
- Geometry modifications: selected main hull; sealed openings; rotated axes and 1:10 scaling; two convex deckhouse envelopes joined with the hull; positioned at the declared tank station; exported contact STL and normals. Small closed reference meshes are included to reproduce the actual inputs. They are not production CAD or a new authoritative vessel model.
- Mass is derived from the visual hull at an assumed draft. KG and radii of gyration are declared assumptions, not measured loading data. Reference dimensions are linked through [NTNU's vessel page](https://www.ntnu.edu/gunnerus/).
- Videos retain the IHB NTNU / vessel.js / MIT credit. Rendering changes camera, lighting, material roughness and water representation; no authorship of the upstream model is claimed.

The repository's GPL code license does not replace the MIT notice for this asset. The Gunnerus permission is independent of, and does not relax, the turtle model's redistribution restriction.
