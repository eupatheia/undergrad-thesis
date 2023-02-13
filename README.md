## Progress Log

- [x] **01.31.23: Color Buddhabrot**
  - color_buddhabrot.c generates the buddhabrot with varying maxIterations for each color channel (also allows variable thread count, regions, size)
- [x] **02.07.23: Mandelbrot Shader**
  - fixed color_buddhabrot.c to use pthread_cond_t instead of pthread_barreir_t (which does not run on macOS)
  - imaged buddhabrots in smaller and smaller regions to confirm self-similarity
  - implemented Mandelbrot visualization with a shader (see mandelbrot.cpp and mandelbrot.fs in the annex repo)
- [x] **02.14.23:**

## Remaining Work Plan

- [ ] Read and implement “Fractal Flames”, 2008, Draves and Reckase
- [ ] Visualizing fractals in three dimensions, e.g.
  - [ ] Mandelbulb
  - [ ] “Strange Attractors: Creating Patterns in Chaos (Chpts 1-3)”.  
  - [ ] Fractal Anatomy: Imaging internal fractal structures (SIGGRAPH paper)
- [ ] Fractals in procedural graphics: clouds and terrain
  - [ ] Texturing and modeling: A Procedural Approach by Ebert, Musgrave, Peachey, Perlin, and Worley (Chpts “A brief introduction to fractals”, “Fractal solid textures: clouds”, and “procedural fractal terrains”
 
Readings: https://www.dropbox.com/sh/ov7c5nmmzrbm3d9/AAAwfBCoOcqTh3RR2El8jKxJa?dl=0
