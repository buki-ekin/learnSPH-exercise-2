# Assignment 5(a) — incremental motion tests

## What was requested

Test zero forces, pressure, viscosity, two-block collision, gravity, floor interactions and boundary viscosity.

## What was implemented

Additional commits: 4fc8040 (gravity), c4b4991 (floor and boundary viscosity). Each stored scene is rendered in full. assignment5_floor is a partial 35-frame diagnostic run, not a completed 3-second run. The dt025 variant tests a smaller maximum time step.

## Exact source revision

[Commit effdd44](https://github.com/buki-ekin/learnSPH-exercise-2/commit/effdd44). Later commits contain cumulative changes; use a detached worktree to build this historical state.

## Videos

- [assignment5_zero](../../showcase/all_scenes/assignment5_zero.mp4)
- [assignment5_pressure](../../showcase/all_scenes/assignment5_pressure.mp4)
- [assignment5_viscosity](../../showcase/all_scenes/assignment5_viscosity.mp4)
- [assignment5_collision](../../showcase/all_scenes/assignment5_collision.mp4)
- [assignment5_freefall](../../showcase/all_scenes/assignment5_freefall.mp4)
- [assignment5_floor](../../showcase/all_scenes/assignment5_floor.mp4)
- [assignment5_floor_1s](../../showcase/all_scenes/assignment5_floor_1s.mp4)
- [assignment5_floor_1s_dt025](../../showcase/all_scenes/assignment5_floor_1s_dt025.mp4)
- [assignment5_floor_friction](../../showcase/all_scenes/assignment5_floor_friction.mp4)

## Raw outputs

Original VTK, PNG and Blender files are available in the [complete output release](https://github.com/buki-ekin/learnSPH-exercise-2/releases/tag/complete-exercise-archive-v1). Folder names are preserved inside each ZIP.
