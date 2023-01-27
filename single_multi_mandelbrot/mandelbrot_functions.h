#ifndef MANDELBROT_FUNCTIONS_H_
#define MANDELBROT_FUNCTIONS_H_

extern void generatePalette(struct ppm_pixel * palette, int maxIterations);
extern void computeMandelbrot(int size, int i, int j, float xmin, float xmax,
    float ymin, float ymax, int maxIterations, struct ppm_pixel ** pixels,
    struct ppm_pixel * palette);

#endif
