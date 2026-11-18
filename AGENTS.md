# Instructions for agents working in this repository

- Read `docs/AGENT_HANDOFF.md` and `docs/REPRODUCE.md` first.
- Keep source and derived turtle geometry private. Never commit model archives, STL/OBJ/FBX/Blender files, sampled hull functions, BEM panels, boundary particle coordinates, roof maps or render caches.
- Keep full solver fields, binaries, toolchains and credentials out of this repository. Preserve the upstream sources, attribution and license notices.
- Preserve accepted reference inputs/results and final videos. New variants get new run names; record changes and actual completion status.
- Potential flow and SPH are independent. The focused-wave intro and regular-wave coasting scene are different conditions. Do not claim validated peak pressures, whole-roof green water, flooding or historical performance.
- The user requested an approval checkpoint before substantial code changes, prolonged debugging and other high-usage implementation work. Routine inspection and documentation can continue; explain a concrete implementation scope before asking approval. Existing explicit authorization takes precedence.
