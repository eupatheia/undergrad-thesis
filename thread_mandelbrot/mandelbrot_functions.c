#include <stdlib.h>
#include <time.h>
#include "read_ppm.h"

// generate palette colors
void generatePalette(struct ppm_pixel * palette, int maxIterations) {
  unsigned int basered, basegreen, baseblue;
  srand(time(0));

  basered = rand() % 255;
  basegreen = rand() % 255;
  baseblue = rand() % 255;
  for (int i = 0; i < maxIterations; i++) {
    palette[i].red = basered + rand() % 100 - 50;
    palette[i].green = basegreen + rand() % 100 - 50;
    palette[i].blue = baseblue + rand() % 100 - 50;
  }
}

// compute mandelbrot set
void computeMandelbrot(int size, int start_row, int end_row, int start_col,
    int end_col, float xmin, float xmax, float ymin, float ymax,
    int maxIterations, struct ppm_pixel ** pixels, struct ppm_pixel * palette) {

  float xfrac, yfrac, x0, y0, x, y, xtmp;
  int iter;

  for (int i = start_row; i < end_row; i++) {
    for (int j = start_col; j < end_col; j++) {
      xfrac = (float) j / size;
      // must flip rows, otherwise -y axis points up
      yfrac = (float) (size - i - 1) / size;
      x0 = xmin + xfrac * (xmax - xmin);
      y0 = ymin + yfrac * (ymax - ymin);

      x = 0;
      y = 0;
      iter = 0;
      while (iter < maxIterations && x*x + y*y < 2*2) {
        xtmp = x*x - y*y + x0;
        y = 2*x*y + y0;
        x = xtmp;
        iter++;
      }
      
      if (iter < maxIterations) {
        // escaped
        pixels[i][j].red = palette[iter].red;
        pixels[i][j].green = palette[iter].green;
        pixels[i][j].blue = palette[iter].blue;
      } else {
        // did not escape, use color black
        pixels[i][j].red = 0;
        pixels[i][j].green = 0;
        pixels[i][j].blue = 0;
      }
    }
  }
}
