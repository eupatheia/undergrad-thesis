# Visualizing Chaotic Systems

This work explores 5 different iterative algorithms that form fractals, all implemented on CPU using threads:
- Filled Julia Set
- Mandelbrot Set
- Mandelbulb
- Buddhabrot Set
- Fractal Flames

## Mandelbulb Demo

The Mandelbulb is a 3D Mandelbrot-like object that has fractal qualities.  Here, we render it using raymarching and signed distance functions (SDFs) to approximate the surface (see [Quilez's Mandelbulb article](https://iquilezles.org/articles/mandelbulb/) for the algorithm).  Parts of the volumetric coloring algorithm were adapted from [this Shadertoy page](https://www.shadertoy.com/view/wt33Wl), and the soft shadow algorithm came from [this article](https://iquilezles.org/articles/rmshadows/).  This implementation outputs a ppm image and does not allow for user interaction.  However, the image size, number of threads, anti-aliasing filter size, and coloring scheme can be adjusted by command line arguments.  Camera position and other parameters currently cannot be adjusted from the command line.  The gif below was compiled from 36 separate images.

![ezgif com-crop-2-2-2-2](https://user-images.githubusercontent.com/75283980/234853201-abf16bc2-806a-49e4-9c38-7fc4289a6919.gif)

### How to Build

Download the `mandelbulb` folder and navigate to the directory in the terminal.

To build the program, type `make mandelbulb`.

To run the program, use the following command, with optional flags:
```
./mandelbulb -s <size> -p <numProcesses> -a <antialiasingFilterWidth> -c <colorSetting>
```
- `size`: the width or height of the output image, in pixels (default 300)
- `numProcesses`: number of threads to use (default 4)
- `antialiasingFilterWidth`: the size of the filter used for anti-aliasing (default 1)
  - e.g. an argument of 2 means the program will average 2x2=4 point samples within a pixel to compute the final color for the pixel
- `colorSetting`: the coloring scheme to be used (default normal)
  -  1 = phong shading
  -  2 = total reflection
  -  3 = total refraction
  -  4 = dielectric
  -  5 = chromatic dispersion
  -  6 = volumetric
  -  (other) = normal
  
### Features and Results

#### Anti-Aliasing

<img src="https://user-images.githubusercontent.com/75283980/234854498-43b57bad-6b69-4c2b-8610-368d24bbdd56.png" height=450px/> <img src="https://user-images.githubusercontent.com/75283980/234854571-b5fef8c3-12ec-43ac-ac80-ce3516df3206.png" height=450px/>

#### Surface Normals as Colors, Camera Rotation

<img src="https://user-images.githubusercontent.com/75283980/234853730-0e4b092d-2f14-43bd-ad4d-6639209fa09f.png" height=300px/> <img src="https://user-images.githubusercontent.com/75283980/234853800-0ee18ce5-7096-4579-9df9-5e558a174bd8.png" height=300px/> <img src="https://user-images.githubusercontent.com/75283980/234853812-31a976ad-f847-4758-b81d-0270fa8d182c.png" height=300px/>

#### Phong Model Shading + Soft Shadows, Gradient Background

<img src="https://user-images.githubusercontent.com/75283980/234855330-3a11d28c-9a90-43b3-ac69-2ef55a9fc72d.png" height=500px/>

#### Total Reflection, Procedural Cube Map

<img src="https://user-images.githubusercontent.com/75283980/234855632-858c6181-18f4-43dd-8023-34f3b1ed5f36.png" height=500px/>

#### Total Refraction

<img src="https://user-images.githubusercontent.com/75283980/234855778-5620971d-c645-4f72-9b86-35d7e30ce6ce.png" height=500px/>

#### Dielectric Material

<img src="https://user-images.githubusercontent.com/75283980/234855845-05cd14be-0175-4523-aa79-0ca0151b9bfa.png" height=500px/>

#### Chromatic Dispersion

<img src="https://user-images.githubusercontent.com/75283980/234856096-14d3ae6a-b88a-4e93-b5b3-a57444296da3.png" height=300px/> <img src="https://user-images.githubusercontent.com/75283980/234856710-f98450d1-58c5-418b-b9a6-bb86d2558bd7.png" height=300px/> <img src="https://user-images.githubusercontent.com/75283980/234856192-047ce936-e2d2-4bbc-912a-6bf446be91c6.png" height=300px/> 

#### Volumetric Rendering

<img src="https://user-images.githubusercontent.com/75283980/234857220-2bd38ea1-e14a-4753-8f16-53b79c3b1d33.png" height=300px/> <img src="https://user-images.githubusercontent.com/75283980/234857238-5e94e11c-64cf-46f3-9799-17c57b2b57c0.png" height=300px/> <img src="https://user-images.githubusercontent.com/75283980/234857253-8d4f2f6f-6df7-40ed-8c69-aed99a3de64f.png" height=300px/>

---

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
- [x] **03.21.23: Fractal Flames IV**
  - iterated another close point (0.000001 more in each direction) to calculate Lyapunov exponent to measure chaos, as defined in Sprott (1993)
    - for two arbitrarily close points, the average log ratio of distances before and after an iteration
    - clamped distances at 0.000001 to prevent undefined division/log errors
  - defined new objective as colorfulness (percent non-black pixels) >= 0.15 and positve Lyapunov exponent
    - "best" specimen is defined first by best colorfulness, then by acceptable Lyapunov exponent
  - only use symmetry 1 across the board (search and random) to get accurate Lyapunov
  - 1+1 mu+lambda still ideal
- [x] **03.28.23: Mandelbulb I**
  - implemented raymarching algorithm on a test scene (see `raymarcher.c`)
    - phong per-pixel shading, 3D SDFs, distance-aided marching (from learning references: [Intro to Raytracing](https://www.scratchapixel.com/lessons/3d-basic-rendering/introduction-to-ray-tracing/implementing-the-raytracing-algorithm.html) and [Raymarching](https://michaelwalczyk.com/blog-ray-marching.html))
- [x] **04.04.23: Mandelbulb II**
  - added shadows to `raymarcher.c` (buggy)
  - initial implementation of Mandelbulb (order 8, 4 iterations) using Hubbard-Douady potential to get fractal SDF in raymarching
    - algorithm and derivations from [Mandelbulb](https://iquilezles.org/articles/mandelbulb/) (Quilez, 2009)
    - finds vague shape centered at origin, still buggy (see `mandelbulb.c`)
- [x] **04.11.23: Mandelbulb III**
  - fixed Mandelbulb (equation and color errors), but only allows camera position along negative z-axis (`mandelbulb.c`)
    - uses phong shading currently
  - fixed shading in `raymarcher.c`, attempted to add soft shadows (buggy)
- [x] **04.18.23: Mandelbulb IV + Debugging**
  - fixed bugs in `flame.c`
    - F3 weight was not being printed
    - Reusing Gaussian noise was incorrect with when called with different (mean, stddev)
    - mutate() was not preserving sum of probabilities = 1
  - with bug fixes, 2+2 mu+lambda search seems ideal, jitter by stddev 0.05, colorfulness cutoff 0.33
  - added 4x4 matrix struct implementation with associated helper functions
  - fixed y-axis being upside down in `julia.c` and `thread_mandelbrot.c`
- [x] **04.27.23: Mandelbulb V**
  - implemented camera rotation (inverse view matrix) in `mandelbulb.c`
  - organized helper functions into different files `ray_functions.c` and `shaders.c`, global to local variables
  - implemented anti-aliasing option (average colors of point samples within a pixel to get final color)
  - allow for 7 different coloring schemes, chosen by command line args
  - implemented environment mapping, i.e. procedural cubemap, gradient background
  - colorings: normal, phong + soft shadows, reflection, refraction, dielectric, chromatic dispersion, volumetric
