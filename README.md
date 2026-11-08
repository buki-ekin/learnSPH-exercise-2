# learnSPH — Exercises 1 and 2

Course framework extended with SPH density estimation, boundary handling, WCSPH forces, simulation experiments and Blender visualization. See the [complete assignment folders](assignments/README.md), [all scene videos](showcase/all_scenes/), and [downloadable raw-output archives](https://github.com/buki-ekin/learnSPH-exercise-2/releases/tag/complete-exercise-archive-v1).

See [Exercise 2 progress](EXERCISE_2.md) for requirements, implementation and the actual commit history.

## Build and run

Run commands from the repository root. Requires a C++17 compiler and CMake; dependencies are included under `extern/`. OpenMP is optional. Blender and FFmpeg are needed for visualization.

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j 4
./build/tests/learnSPH_tests --reporter compact
./build/app/learnSPH_density
./build/app/learnSPH_app
./build/app/learnSPH_custom
```

- `learnSPH_density`: static density, boundary and force experiments.
- `learnSPH_app`: Assignment 5(b), dam-break reference snapshots, 0–0.38 s at 50 FPS; output `output/assignment5_dambreak_reference`.
- `learnSPH_custom`: Assignment 5(c), two opposing blocks, 24,000 initial particles, 3 s at 50 FPS; output `output/assignment5_custom`.
- `learnSPH_benchmark` and `learnSPH_experiments`: earlier neighbor-search and integration experiments.

## Videos and rendering

Ready-to-view videos are in [showcase](showcase/). Original VTK files, individual PNGs and Blender snapshots are preserved as ZIP assets in the complete-output release. Build products remain excluded from Git.

On macOS:

```bash
/Applications/Blender.app/Contents/MacOS/Blender --background --python blender/render_custom.py -- output/assignment5_custom output/assignment5_custom_blender
ffmpeg -framerate 15 -start_number 0 -i output/assignment5_custom_blender/custom_%04d.png -frames:v 151 -c:v libx264 -crf 18 -pix_fmt yuv420p -movflags +faststart output/assignment5_custom_blender/custom_scene_3d.mp4
```

On other platforms replace the Blender executable path. The custom renderer handles changing particle counts and uses a shared speed color scale. Playback is 15 FPS: 151 frames give approximately 10.07 seconds, showing 3 seconds of simulation in slow motion. `custom_final.blend` stores only the final static scene.

## Source layout

- `learnSPH/`: kernel, sampling, neighbor search, integration and VTK I/O.
- `app/density.cpp`: Assignments 1–4 experiments.
- `app/main.cpp`: final Assignment 5(b) configuration.
- `app/custom_scene.cpp`: Assignment 5(c).
- `blender/`: visualization scripts.
- `tests/`: automated tests inherited and extended during the exercises.

Earlier Assignment 5(a) configurations are preserved in Git history rather than exposed as runtime switches. Existing comments may still mention earlier experiments; the executable configuration is authoritative. The optional moving-boundaries assignment is not implemented. Third-party components retain their own licenses in `extern/` and `tests/`.
