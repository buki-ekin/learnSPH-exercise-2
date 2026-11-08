# Exercise 2 — requirements and implementation history

These are the real incremental commits, not a reconstructed history. Assignment 1's initial commit also includes the existing Exercise 1 framework.

| Assignment | Requested work | Implementation and evidence | Commit |
|---|---|---|---|
| 1(a–b) | Estimate SPH density; sample a one-litre cube with 20³ particles and examine the distribution | Constant mass, self contribution and neighbor sum; 8,000 particles, total mass 1 kg; density histogram | `7a77d86` |
| 2(a–b) | Sample boundaries and compute boundary volumes/masses | Triangle sampling; inverse kernel-sum volumes including self; mass from rest density | `5903e51` |
| 3(a–b) | Include boundary density and visualize its effect | Fluid plus boundary density; VTK scalar output and Blender comparison | `8647c7f` |
| 4(a) | Implement equations of motion | Clamped linear pressure, fluid/boundary pressure and viscosity, gravity; static checks | `35cb85f` |
| 5(a), baseline–iii | Incremental zero-force, pressure, viscosity and two-block collision tests | Baseline tests and colliding blocks; Blender snapshots | `effdd44` |
| 5(a), iv | Single block under gravity | Free-fall check | `4fc8040` |
| 5(a), v–vi | Floor interaction and boundary viscosity | Floor pressure and damping comparison | `c4b4991` |
| 5(b) | Reproduce the dam-break case and color by velocity | Nominal particle mass; specified fluid dimensions and coefficients; reference snapshots and speed coloring | `1bcbbca` |
| 5(c) | Own scene, ≥20,000 particles, several seconds, ≥10 s video | Two opposing blocks, 24,000 particles, 3 s simulation, particle deletion, Blender rendering | `e14ca3a` |

## Numerical method

Cubic spline kernel, h = 1.2d, support radius 2h, particle diameter d = 2r. Fluid mass is rho0*d³. Boundary volumes are computed from boundary kernel sums. Pressure uses max(0, B*(rho-rho0)), as in the course's simplified model. Accelerations combine fluid and boundary pressure, viscosity and gravity; integration is semi-implicit Euler with bounded CFL time steps.

For 5(b): fluid dimensions 0.15 × 0.25 × 0.50 m; radius 0.005 m; rest density 1000 kg/m³; stiffness 1000; fluid viscosity 0.0025; boundary viscosity 0. The current maximum step is 0.05 ms, below the assignment's 0.25 ms ceiling. The tank implementation adds one particle diameter on each side of its nominal dimensions. This is an interpretation of the assignment's enlargement wording, not a verified unique interpretation.

## Custom scene and observed results

Two blocks of 12,000 particles each start with y velocities +0.5 and -0.5 m/s in an open-top tank. The recorded 3-second run ended with 23,814 particles: 186 (0.775%) were removed. Final-frame maximum speed was 0.227049 m/s, and maximum density approximately 1005.22 kg/m³. Final diagnostics reported no particles above the rim or outside the tank footprint. These observations are not a proof of stability for arbitrary settings.

The main challenges were keeping persistent particle arrays aligned after deletion, updating CompactNSearch, and rendering snapshots with changing particle counts. The renderer rebuilds the point mesh per snapshot and preserves a shared speed scale and fixed oblique camera. The final maximum-speed and displacement summaries describe surviving particles in the final state, not global maxima over the whole run.

## Revisit an earlier step

Use `git show <commit>` to inspect a step. To build an earlier state without altering current work:

```bash
git worktree add --detach ../learnSPH-step effdd44
```

Build that directory separately. Later source files contain cumulative changes. The optional moving-boundaries extension has not been implemented.
