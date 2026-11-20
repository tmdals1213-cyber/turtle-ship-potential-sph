# Instructions for agents working in this repository

- Read `docs/AGENT_HANDOFF.md` and `docs/REPRODUCE.md` first.
- For case 06, read `docs/TANKER_REPRODUCE.md` and `docs/TANKER_MODEL_NOTES.md`. Its analytic tank/truck geometry is generated locally; no restricted ship assets are used. Keep immutable compact references and public final films, not full fields. Native Chrono is required; the older ship build has it disabled.
- Keep source and derived turtle geometry private. Never commit model archives, STL/OBJ/FBX/Blender files, sampled hull functions, BEM panels, boundary particle coordinates, roof maps or render caches.
- The explicit exception is the allowlisted MIT-derived **Gunnerus** reference contact geometry under case 05, with its upstream license notice. This does not permit publishing any turtle geometry, generated particle fields or render caches. Read `docs/GUNNERUS_REPRODUCE.md` and `docs/MARINE_MODEL_NOTES.md` for this case; `gunnerus/` is its separate workflow.
- Keep full solver fields, binaries, toolchains and credentials out of this repository. Preserve the upstream sources, attribution and license notices.
- Preserve accepted reference inputs/results and final videos. New variants get new run names; record changes and actual completion status.
- Potential flow and SPH are independent. The focused-wave intro and regular-wave coasting scene are different conditions. Do not claim validated peak pressures, whole-roof green water, flooding or historical performance.
- The user requested an approval checkpoint before substantial code changes, prolonged debugging and other high-usage implementation work. Routine inspection and documentation can continue; explain a concrete implementation scope before asking approval. Existing explicit authorization takes precedence.
