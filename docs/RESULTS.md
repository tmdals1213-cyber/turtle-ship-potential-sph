# Results and limits

This package preserves the accepted visualization demonstration and compact numerical records. It is not experimental validation, a historical performance prediction, or a solver ranking.

## Physical assumptions

The scenario uses 1:10 Froude scaling, water density 1,025 kg/m3 and a 154.4 kg model mass (154.4 t nominal full-scale equivalent). The mass budget combines a published 140.4 t reconstruction self-weight with illustrative crew, equipment and stores. The 26.06 m reference hull length is derived by subtracting a reported 9.21 m tail from 35.27 m overall length; it is not a published Lpp or Lwl. CG, inertia and loading positions are assumptions. Froude similarity does not also provide full-scale Reynolds similarity.

Water acts on a closed, simplified hull that includes an initially dry upper boundary. Decorative mast, oars and rigging are display geometry. Air, wind, trapped air, interior flooding and structural damage are absent. Public physical properties are in `study/results/turtle_history_model.json`; source geometry remains private.

## Potential-flow result

The 6,165-panel BEM result has heave amplitude 0.1028 m and pitch amplitude 12.55 degrees for the selected wave. The large pitch is outside a comfortable small-motion regime, so this result is qualitative. Panels remain on the mean wetted hull; a moving display waterline does not mean the hydrodynamic mesh was remeshed every frame. CPU response solve: 153.85 s; including surface evaluation: 213.79 s.

## Coasting SPH result

The 18 x 2.4 m tank has 1 m water depth. The initial body speed is +0.5 m/s with no continuing thrust. Regular-wave input H=0.24 m, wavelength=3 m; measured fundamental wave heights were about 0.176–0.188 m without incident/reflected separation. The 943,816-particle calculation completed 16.0001 model seconds in 688.30 s on RTX 5080, with zero excluded particles. Last saved motion time: 15.9001 s. Pitch: about -6.35 to +7.32 degrees. The early 2–3 cm mean-heave bias remains a discretization/startup limitation.

## Focused-wave result

The tank is 24 x 4.8 m with 1.2 m water depth. The hull starts at rest. First-order JONSWAP/NewWave piston components use a 0.65 m linear coherent crest target, 2.2 s peak period, focus x=1 m and target time 7.5 s. These are input parameters: the nonlinear crest developed earlier. The empty flume measured a maximum elevation near 0.568 m at x=1 m around t=4.99 s.

The 2,885,421-particle ship case completed 10.0002 model seconds in 1,552.27 s, 42,028 native steps and zero excluded particles. Heave range: -0.181 to +0.385 m; relative Euler-Y pitch: -18.04 to +21.29 degrees. Pitch signs follow the stored rotation convention, not an unstated naval sign convention. Intro interval: t=4.60–5.92 s, rendered into three video seconds.

The central-roof probe detected no water above its footprint. A separate upper-hull probe detected a maximum of 102 particles, mainly near the forward upper surface. This supports a bow-impact/upper-hull-contact scene, **not** a claim of a water sheet covering the entire roof.

One late isolated parcel reached 11.41 m/s at about t=8.40 s; only two parcels exceeded 6 m/s in that snapshot, and the final snapshot had none above 6 m/s. The intro interval's stored maximum was about 5.32 m/s. The density range over saved snapshots was 972.49–1,068.73 kg/m3. See the included `late_peak_audit.json`; neither isolated peak velocity nor impact pressure has been independently validated. Zero excluded particles alone is not proof of accuracy.

## Visualization

Rigid poses come from computed floating particles. The intro follows the same water-particle IDs between saved outputs, then reconstructs a smoothed 3-D occupancy surface. Dilute displayed parcels are sampled simulated water, not invented foam. Small droplets and pressure peaks are not resolved claims. The displayed water region is a crop of the computational tank. The potential and SPH segments solve different mathematical problems on different hardware; their timings are not a speedup ratio.
