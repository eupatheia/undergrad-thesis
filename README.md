## Progress Log

- [x] **01.31.23: Color Buddhabrot**
  - `color_buddhabrot.c` generates the buddhabrot with varying `maxIterations` for each color channel (also allows variable thread count, regions, size)
- [x] **02.07.23: Mandelbrot Shader**
  - fixed `color_buddhabrot.c` to use `pthread_cond_t` instead of `pthread_barrier_t` (which does not run on macOS)
  - imaged buddhabrots in smaller and smaller regions to confirm self-similarity
  - implemented Mandelbrot visualization with a shader (see `mandelbrot.cpp` and `mandelbrot.fs` in the annex repo)
- [x] **02.14.23: Julia Sets and Attractors I**
  - drew one-parameter complex quadratic Julia sets as ppm image (see `julia.c`) and with shader (see `julia.cpp` and `julia.fs` in annex repo)
    - coloring based on how quickly points diverge
  - added same coloring scheme to Mandelbrot shader
  - read Ch. 1-4 of [*Strange Attractors*](https://sprott.physics.wisc.edu/SA.HTM) (Sprott 1993)
- [x] **02.21.23: Fractal Flames I**
  - read ["Fractal Flames"](https://flam3.com/flame_draves.pdf) (Draves and Reckase, 2008)
  - implemented spherical two-function system with log density display, attempt to add color (see `flame.c`)
- [x] **02.28.23: Fractal Flames II**
  - fixed coloring algorithm
  - added density estimation supersampling with Gaussian kernels
  - used multiple threads for rendering step (coloring pixels using counts)
  - approximated example flames (see `fractal_flames/spherical` and `fractal_flames/linear-swirl-spiral`, cf. Figures 3 and 5 in Draves 2008)
  - restructured code to allow for random generation of the iterative function system:
    - number of functions (2-4)
    - symmetry (1-3)
    - choice of functions (among 17 possible)
    - weights of functions (i.e. probability of being chosen while iterating)
    - affine transform coefficients (all in range [-1,1])
    - color palette (among 12 possible)
- [x] **03.14.23: Fractal Flames III**
  - read Ch. 2 of *Procedural Content Generation in Games*, ["The Search-based Approach"](https://link.springer.com/content/pdf/10.1007/978-3-319-42716-4_2.pdf)  (Togelius and Shaker, 2016)
  - implemented mu+lambda evolutionary search on fractal flames using percentage of non-black pixels (nonzero counts) and 0.3 threshold for "interest"
    - only symmetry 1 used in search to prevent objective function bias towards high symmetry to fill space
    - search usually stops after 1 (or 2) generations with mu = 1 and lambda = 1, not much improvement over random generation
  - restructured code to handle random and search options (see `flame.c`)
  
## Remaining Work Plan

- [ ] Visualizing fractals in three dimensions, e.g.
  - [ ] Mandelbulb
  - [ ] “Strange Attractors: Creating Patterns in Chaos (Chpts 1-3)”.  
  - [ ] Fractal Anatomy: Imaging internal fractal structures (SIGGRAPH paper)
- [ ] Fractals in procedural graphics: clouds and terrain
  - [ ] Texturing and modeling: A Procedural Approach by Ebert, Musgrave, Peachey, Perlin, and Worley (Chpts “A brief introduction to fractals”, “Fractal solid textures: clouds”, and “procedural fractal terrains”
- [ ] Buddhabrot equivalent for non-quadratic functions?
 
Readings: https://www.dropbox.com/sh/ov7c5nmmzrbm3d9/AAAwfBCoOcqTh3RR2El8jKxJa?dl=0
