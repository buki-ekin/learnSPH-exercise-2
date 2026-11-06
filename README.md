# learnSPH - Exercise 1

This project completes Exercise Sheet 1 with a simple particle cube. It does not
use Orion, Gemini or another mesh. The whole cube volume is filled with particles.

## Build and test

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j 6
./build/tests/learnSPH_tests --reporter compact
```

## Make the VTK frames

```bash
cd build
./app/learnSPH_app
cd ..
```

The program creates 91 files in `build/output/vtk`. The animation has 30 frames
per second. Every file stores position, height and velocity data.

## Run the experiments

```bash
./build/app/learnSPH_benchmark
./build/app/learnSPH_experiments
```

The benchmark compares the simple search with CompactNSearch. The experiment
compares Euler with the exact result and shows tunneling with and without CFL.

## Render with Blender

```bash
/Applications/Blender.app/Contents/MacOS/Blender \
  --background \
  --python blender/render_cube.py \
  -- build/output/vtk build/output/blender

ffmpeg -y -framerate 30 \
  -i build/output/blender/frame_%04d.png \
  -c:v libx264 -pix_fmt yuv420p \
  build/output/blender/cube_particles.mp4
```

The Blender script makes small spheres, adds a floor, camera and lights, and
colors the particles blue, purple, pink, orange and yellow using their position.

## Important files

- `learnSPH/kernel.*`: cubic spline kernel and gradient
- `learnSPH/neighborhood.*`: simple neighbor search and particle grid
- `learnSPH/simulation.*`: Euler, CFL and AABB removal
- `learnSPH/io.*`: VTK writing
- `tests/sample_tests.cpp`: automatic tests
- `app/main.cpp`: falling particle cube
- `app/benchmark.cpp`: neighbor-search timing
- `app/experiments.cpp`: Euler and tunneling examples
- `blender/render_cube.py`: PNG rendering
