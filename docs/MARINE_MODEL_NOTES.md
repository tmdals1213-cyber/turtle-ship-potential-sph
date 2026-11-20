# Marine SPH: what the model represents

These notes distinguish the **completed Gunnerus demonstration** from possible extensions. It predicts the motion of the declared coarse model; it is not a validated prediction of the real vessel's resistance, capsize limit or impact pressure. [Reference setup and results](../SPH_CFD-incomp-free-surface-6DOF_GPU_Cpp-CUDA-Python_05-gunnerus-focused-wave).

## A ship is one rigid body, not independent solid particles

The solid boundary samples locate fluid contact. Their forces and moments are accumulated for a single body with total mass, centre of gravity (CG) and inertia tensor. Its translation and rotation move the boundary samples together. This differs from a DEM assembly whose grains can rearrange through contacts. See the [DualSPHysics rigid-body formulation](https://github.com/DualSPHysics/DualSPHysics/wiki/3.-SPH-formulation#38-fluid-driven-objects).

For a rigid vessel with fixed internal contents, steel, engines and cargo can be represented through **M, CG and I**. Filling the entire external hull volume with steel density would give the wrong mass. The wetted external shape determines the pressure/buoyancy distribution; CG determines its moment arm, and inertia determines angular acceleration. Equal hull shape and mass do not imply equal motion when CG or inertia changes.

| Current input | Model value / basis |
|---|---|
| Scale | 1:10, gravity 9.81 m/s², reference water density 1,025 kg/m³ |
| Mass | 398.8999 kg; closed visual-hull displacement at an assumed 0.27 m draft |
| Full-scale equivalent mass | 398.9 t; **not a measured vessel displacement** |
| Initial CG, tank coordinates | (5.3497144, approximately 0, 0.13) m; z=0 is still water |
| KG | Assumed 4.0 m full equivalent; longitudinal CG initially aligned with buoyancy centre |
| Principal inertia about CG | (45.1390, 286.5123, 286.5123) kg·m²; assumed radii of gyration |
| Six degrees of freedom | All free; initial +x speed and an assumed world-fixed +x force |
| Contact surface | Sealed hull plus two simplified deckhouse envelopes; masts, crane and rails are visual only |

Adding the deckhouse contact envelopes did not add a second structural mass. Flooding, sloshing, moving cargo, flexible hull response and air loads require additional models. The water-only contact geometry can become submerged above its initial waterline. The full decorative visual mesh is not the contact surface.

## Why Froude similarity is used

Gravity and inertia dominate wave propagation and the large vessel motions targeted here. Preserve `Fr = U / sqrt(g L)` when changing scale. With length ratio λ = full/model = 10, equal gravity and reference density give:

| Quantity | Model → full-scale equivalent |
|---|---|
| Length, wave height, depth | ×10 |
| Speed, time, wave period | ×√10 |
| Mass, force | ×1,000 |
| Moment | ×10,000 |
| Mass moment of inertia | ×100,000 |
| Pressure | ×10 |
| Angle, linear acceleration | unchanged |

Here 1.30145 m/s corresponds to 8 kn, the 2.6 s component peak period corresponds to 8.22 s, and 16 model seconds correspond to 50.6 equivalent seconds. The film displays **model time and model-unit loads**. Initial Fr based on Lpp=3.39 m is approximately 0.226. These are similarity conversions, not validation of the real ship. See [ITTC resistance-test definitions and scaling](https://ittc.info/media/1816/75-02-02-01.pdf).

Reynolds similarity is not simultaneously preserved with the same physical kinematic viscosity: `Re_model / Re_full = λ^(-3/2) = 1/31.62`. Matching both would require a different model viscosity. This run instead uses SPH artificial viscosity α=0.01; it is **not** a wall-resolved water-friction calculation. Full-scale friction/resistance cannot be recovered simply by multiplying the current force by 1,000.

## What SPH, FVM, IBM and DEM each mean

| Approach | Relevant strength | Requirement / limitation for this task |
|---|---|---|
| Body-fitted FVM, usually RANS + VOF | Established ship resistance and boundary-layer workflows | Moving/overset meshes and interface treatment need suitable resolution; breaking-wave problems remain possible |
| IBM on a background grid | Moving solid geometry without fitting/deforming the whole grid | Solid-boundary accuracy and a free-surface model are still needed |
| SPH | Moving fluid samples naturally follow large free-surface deformation and fragmentation | Kernel support, dissipation, wall treatment, spacing and time-step sensitivity still need verification |
| DEM | Contact and rearrangement of discrete solid grains or bodies | Not a substitute for the continuum water solver |

FVM is a discretization, IBM a boundary treatment, and VOF an interface method: **FVM + IBM + VOF can coexist**. The comparison slide is illustrative, not three matched CFD runs or a performance benchmark. SPH particles represent fluid volumes, not physical water molecules/droplets. See [fixed-grid immersed-boundary formulation](https://engineering.jhu.edu/fsag/wp-content/uploads/2014/06/jcp-vol174.pdf), [OpenFOAM solver families](https://www.openfoam.com/documentation/user-guide/a-reference/a.1-standard-solvers), and [DualSPHysics formulation](https://github.com/DualSPHysics/DualSPHysics/wiki/3.-SPH-formulation).

SPH is useful for breaking-wave interaction and large vessel motion, but the best method depends on the quantity being measured. For accurate wall friction, propeller thrust/torque or efficiency, a verified RANS workflow is often the practical baseline. Published [SPH planing-hull resistance](https://www.mdpi.com/2077-1312/9/3/244) and [seakeeping studies](https://www.mdpi.com/2077-1312/11/4/700) show that quantitative SPH work is possible with appropriate verification; they do not validate this hull or setup.

## Particle count and time resolution

The completed run used 9,985,001 total particles: 8,581,760 fluid, 13,638 floating boundary, 119,400 moving piston and 1,270,203 fixed boundary. Uniform initial spacing is 0.05 m model scale (0.5 m full equivalent): only about **19.2 intervals across the beam and 5.4 through the draft**. A large total count does not imply a finely resolved hull boundary layer or thin water film.

Adaptive solver steps were 0.0789–0.2551 ms, mean approximately 0.220 ms, over approximately 72,750 steps. Particle fields were saved every 0.1 model seconds; native body monitors every 0.005 s. Video frames use temporal interpolation and are not extra solver steps. WCSPH time stepping is limited by acoustic/velocity, acceleration and dissipative constraints; see [upstream time integration](https://github.com/DualSPHysics/DualSPHysics/wiki/3.-SPH-formulation).

Halving uniform spacing in 3-D implies roughly eight times as many particles and, under similar acoustic constraints, roughly twice as many steps: an order-of-magnitude planning estimate of 16× work, not a measured runtime guarantee. The final strong-wave case has **no systematic space/time independence study**. Fine spray, drainage ports and centimetre-scale sheets are unresolved. Finite fields and zero lost particles alone do not establish accuracy: sampled density ranged from about 910 to 1,144 kg/m³ around the 1,025 reference.

## Refinement near the ship is possible, but not enabled here

DualSPHysics 5.4 includes GPU variable resolution and tracking of moving objects. It uses coupled resolution zones and buffer reconstruction. The checked driver synchronizes zones using the minimum predicted step; this is not local time subcycling. See [upstream changes](https://github.com/DualSPHysics/DualSPHysics/blob/ef3721a861fda961f0e2f9ec4cd317b19de99086/CHANGES.txt), [VRes XML format](https://github.com/DualSPHysics/DualSPHysics/blob/ef3721a861fda961f0e2f9ec4cd317b19de99086/doc/xml_format/_FmtXML__VRes.xml), and [time-step synchronization](https://github.com/DualSPHysics/DualSPHysics/blob/ef3721a861fda961f0e2f9ec4cd317b19de99086/src/source/JSphVResDriver.h).

A possible next study would retain 0.05 m far-field spacing and test a 0.025 m body-tracking region while checking incident-wave preservation, loads and motions. **Not implemented or run for this case.** The current renderer assumes uniform spacing and persistent IDs; variable particle volumes and changing IDs require corresponding changes.

## RANS/LES performance maps can provide SPH inputs

An offline workflow can tabulate resistance `R(U)`, propeller thrust `T(Va,n)` and torque `Q(Va,n)`, then query those maps using the evolving SPH state. RANS is a reasonable starting point for mean performance; LES is an option when resolved unsteadiness is the question. This is a proposed extension, not the current implementation. [ITTC ship-CFD guidance](https://ittc.info/media/11958/75-03-02-03.pdf) and a [body-force propeller formulation](https://www.mdpi.com/2076-3417/12/21/11273) provide background.

Avoid double-counting resistance already represented by SPH. Add only missing components, or validate a correction against matched calm-water reference cases. Transfer all inputs consistently between model and full scale. An actuator disk applies momentum to the water and an equal/opposite reaction to the ship at the propulsor position and orientation; an external CG force alone creates no propeller wake. Emergence, ventilation and oblique inflow need suitable corrections or additional models.

**Current run:** initial 8 kn equivalent speed, plus a world-fixed 30 N model force at CG ramped from 0.3–1.3 s. It is assumed, not predicted from RPM. No actuator disk, propeller torque, wake, RANS/LES coupling or speed controller is implemented. The ship can be driven backward by the wave.

## Motion is allowed; capsize is not demonstrated

No prescribed roll/pitch cap limits were imposed, but no capsize occurred in the completed 16 s: maximum tilt was about 23.3°. Extrapolating the video is not a capsize prediction. Large coupled rotations should be separately verified for the selected native rigid-body integrator (`RigidAlgorithm=0`); angular integration and orientation handling need benchmark checks before relying on an overturning result. Internal flooding and cargo shifts are absent.

The side/front colors are water pressure sampled 0.025 m outside the contact faces, with unsupported samples shown gray. They are not validated wall traction, pressure-integrated total loads or resolved slam peaks. The native integrated rigid-body forces and moments are reported separately; the pitch moment is not a flexible-hull bending moment.
