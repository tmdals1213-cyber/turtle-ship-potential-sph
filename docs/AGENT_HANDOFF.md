# Handoff to another local agent

Start with README, `cases.json`, `docs/REPRODUCE.md`, `docs/RESULTS.md` and THIRD_PARTY_NOTICES. The accepted deliverable is the 35 s portrait film under root `media/`. Its first 3 s use a separate focused-wave run; the remainder presents potential flow then free-coasting SPH. English captions are intentionally short.

The public package contains Python orchestration, exact reference case inputs and compact results. DualSPHysics C++/CUDA and Capytaine are external dependencies, pinned/documented in THIRD_PARTY_NOTICES. Numerical solver kernels were not modified. BEM and SPH were not coupled.

## Local setup

1. Read AGENTS.md and the current user's instructions.
2. Install the external tools and isolated Python dependencies. Set `DUALSPHYSICS_BIN`; set `TURTLE_ASSET_DIR` if the assets already exist on the receiving machine.
3. Obtain the turtle source from its creator if absent. Do not fetch it from an unofficial mirror or expect it in this archive. Reconstruct local assets using the documented pipeline.
4. Start with no-asset empty-flume input generation. Inspect dimensions, piston motion, particle exclusions, wave gauges and solver completion before a large run.
5. Preserve top-level case snapshots. Put new runs under `study/sph/`, with distinct case names. Do not replace a published result with a failed or incomplete run.

`study/results/turtle_history_model.json` is a numerical seed. Geometry preparation recomputes it; top-level `01-potential/setup/physical_model.json` preserves the original. XML geometry paths in reference folders are placeholders; use the generators, not unresolved XML. Render caches and full fields must be regenerated locally. No credentials or machine-specific user paths are provided.

## Pending work / approval boundary

The accepted turtle demonstration is complete. Candidate follow-on work is a more industrial ship propulsion/resistance example, then possibly thermal/phase-change SPH. Neither is represented here as a completed simulation. Research and light result organization may continue; ask the user before substantial code changes, prolonged debugging or other high-usage implementation work, as requested on 2026-11-18. Describe the concrete scope and expected effort first. Do not launch a new expensive GPU study merely because a repository was cloned.

Do not tune a result to look favorable. Preserve the coarse hydrostatic bias, finite-tank reflections and late isolated velocity observation in the record. Before claiming quantitative resistance, green-water loads or historical survivability, new verification/validation work would be required.
