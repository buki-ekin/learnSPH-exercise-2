# All stored dynamic scenes

Each video uses every stored VTK snapshot from its source folder. The two previously rendered final videos are retained. The dam-break reference movie repeats frames (20 source snapshots, 120 encoded frames). Short diagnostic trials are slowed to four seconds for inspection; no intermediate simulation states are synthesized. Other new videos play at 15 FPS. Color maps actual speed with a fixed range within each scene, not a shared range across different scenes. Cameras remain fixed within each scene.

Static zero-force, pressure and viscosity baselines may look motionless: their stored snapshots are motionless tests. `assignment5_floor` is a partial diagnostic run. Dam-break diagnostic variants are retained for completeness, not presented as validated final configurations.

| Scene | Stored snapshots | Simulated interval (s) | Video duration (s) |
|---|---:|---:|---:|
| [assignment5_collision](assignment5_collision.mp4) | 91 | 3.00 | 6.07 |
| [assignment5_custom](assignment5_custom.mp4) | 151 | 3.00 | 10.07 |
| [assignment5_dambreak](assignment5_dambreak.mp4) | 31 | 1.00 | 4.00 |
| [assignment5_dambreak_dt005](assignment5_dambreak_dt005.mp4) | 13 | 0.40 | 4.00 |
| [assignment5_dambreak_nominal_mass](assignment5_dambreak_nominal_mass.mp4) | 13 | 0.40 | 4.00 |
| [assignment5_dambreak_reference](assignment5_dambreak_reference.mp4) | 20 | 0.38 | 4.00 |
| [assignment5_floor](assignment5_floor.mp4) | 35 | 1.13 | 4.00 |
| [assignment5_floor_1s](assignment5_floor_1s.mp4) | 31 | 1.00 | 4.00 |
| [assignment5_floor_1s_dt025](assignment5_floor_1s_dt025.mp4) | 31 | 1.00 | 4.00 |
| [assignment5_floor_friction](assignment5_floor_friction.mp4) | 31 | 1.00 | 4.00 |
| [assignment5_freefall](assignment5_freefall.mp4) | 91 | 3.00 | 6.07 |
| [assignment5_pressure](assignment5_pressure.mp4) | 91 | 3.00 | 6.07 |
| [assignment5_viscosity](assignment5_viscosity.mp4) | 91 | 3.00 | 6.07 |
| [assignment5_zero](assignment5_zero.mp4) | 91 | 3.00 | 6.07 |

## Re-render a scene

From the repository root (use your Blender executable path):

```bash
blender --background --python blender/render_archive_scene.py -- output/assignment5_collision output/archive_collision
ffmpeg -framerate 15 -i output/archive_collision/custom_%04d.png -c:v libx264 -crf 18 -pix_fmt yuv420p -movflags +faststart collision.mp4
```

Original raw data: [download index](../../assignments/RAW_OUTPUTS.md). Static density comparison and Exercise 1 video: [assignment index](../../assignments/README.md).
