# Publication checks — 2026-11-18

- Parsed every included Python module without syntax errors; parsed included JSON and XML.
- Checked the allowlisted package for restricted geometry, solver binaries, field archives, machine-specific user paths and common GitHub credential markers.
- Verified both MP4 SHA-256 values against their accepted sidecars. The final film is the unchanged 1080 x 1920, 30 fps, 35 s deliverable; its recorded full decode has 1,050 frames and no errors.
- Executed the portable empty-flume input generator without turtle assets. This is an input-generation smoke check, not a new CFD validation.
- Full turtle CFD and film rendering were not repeated for publication; they require external solver installations, private assets and regenerated raw fields.

`MANIFEST.sha256` covers all published files except itself. The transfer archive is made from the committed Git tree, so ignored local check outputs are excluded.
