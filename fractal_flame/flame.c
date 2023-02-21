#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <assert.h>
#include <time.h>
#include <pthread.h>
#include <sys/time.h>
#include <math.h>
#include "functions.h"
#include "read_write.h"

// Iterate through the function system and accumulate counts
void iterate(int size, float xmin, float xmax, float ymin, float ymax, 
    struct point p, float r, float g, float b, int iterations,
    struct pix_counts ** counts, struct ppm_pixel * palette, int * maxCount) {

  // associated color for F_0 is at index 1
  struct ppm_pixel c0 = palette[1];
  float r0 = c0.red / 255.0;
  float g0 = c0.green / 255.0;
  float b0 = c0.blue / 255.0;
  // associated color for F_1 is at index 0
  struct ppm_pixel c1 = palette[0];
  float r1 = c1.red / 255.0;
  float g1 = c1.green / 255.0;
  float b1 = c1.blue / 255.0;

  int k, yrow, xcol;
  for (int i = 0; i < iterations; i++) {
    k = rand() % 2;  // both have probability 0.5, here
    if (k == 0) {
      p = affine(p, 0.562482, -0.539599, -0.42992, 0.397861, 0.501088, -0.112404);
      r = (r + r0) / 2;
      g = (g + g0) / 2;
      b = (b + b0) / 2;
    } else {  // k == 1
      p = affine(p, 0.830039, 0.16248, 0.91022, -0.496174, 0.75046, 0.288389);
      r = (r + r1) / 2;
      g = (g + g1) / 2;
      b = (b + b1) / 2;
    }
    p = spherical(p);
    // do not plot first 20 iterations
    if (i >= 20) {
      // calculate row and col of this point
      // note that this will write the +x, +y quadrant in the bottom right
      // i.e. y-axis is flipped from regular Cartesian orientation
      yrow = round(size * (p.y - ymin) / (ymax - ymin));
      xcol = round(size * (p.x - xmin) / (xmax - xmin));
      if (yrow < 0 || yrow >= size || xcol < 0 || xcol >= size) {
        continue; // out of range
      }
      // increment counters
      counts[yrow][xcol].countR += r;
      counts[yrow][xcol].countG += g;
      counts[yrow][xcol].countB += b;
      counts[yrow][xcol].countA++;
      // update max
      if (counts[yrow][xcol].countA > *maxCount) {
        *maxCount = counts[yrow][xcol].countA;
      }
    }
  }
}

void computeColorBinary(struct ppm_pixel ** pixels, struct pix_counts ** counts,
    int * maxCount, int size) {
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      if (counts[i][j].countA > 0) {
        pixels[i][j].red = 255;
        pixels[i][j].green = 255;
        pixels[i][j].blue = 255;
      } else {
        pixels[i][j].red = 0;
        pixels[i][j].green = 0;
        pixels[i][j].blue = 0;
      }
    }
  }
}

void computeColorLog(struct ppm_pixel ** pixels, struct pix_counts ** counts,
    int * maxCount, int size) {
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      int intensity = (log(counts[i][j].countA) / log(*maxCount)) * 255;
      pixels[i][j].red = intensity;
      pixels[i][j].green = intensity;
      pixels[i][j].blue = intensity;
    }
  }
}

void computeColor(struct ppm_pixel ** pixels, struct pix_counts ** counts,
    struct ppm_pixel * palette, int * maxCount, int size) {
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      float factor = log(counts[i][j].countA) / counts[i][j].countA;
      pixels[i][j].red = palette[(int) (counts[i][j].countR * factor) * 255].red;
      pixels[i][j].green = palette[(int) (counts[i][j].countG * factor) * 255].green;
      pixels[i][j].blue = palette[(int) (counts[i][j].countB * factor) * 255].blue;
    }
  }
}

// helper function to set all struct variables to zero
// for ease of incrementation
void initializeCounts(struct pix_counts ** counts, int size) {
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      counts[i][j].countR = 0.0;
      counts[i][j].countG = 0.0;
      counts[i][j].countB = 0.0;
      counts[i][j].countA = 0;
    }
  }
}

//===========================================================//
//===========================================================//
//===========================================================//

int main(int argc, char* argv[]) {
  int size = 900;
  // use bi-unit square
  float xmin = -1.0;
  float xmax = 1.0;
  float ymin = -1.0;
  float ymax = 1.0;
  int iterations = 9200000;
  struct ppm_pixel ** pixels = NULL;
  struct pix_counts ** counts = NULL;
  struct ppm_pixel * palette = NULL;
  double timer;
  struct timeval tstart, tend;
  char new_file[100];
  float x, y;
  int k;  // to randomly choose a function
  int maxCount = 0;

  int opt;
  while ((opt = getopt(argc, argv, ":s:")) != -1) {
    switch (opt) {
      case 's': size = atoi(optarg); break;
      case '?': printf("usage: %s -s <size>\n", argv[0]); break;
    }
  }
  printf("Generating flame with size %dx%d\n", size, size);

  // allocate memory for pixels
  pixels = malloc(sizeof(struct ppm_pixel *) * size);
  if (pixels == NULL) {
    printf("Error: failed malloc.  Exiting...\n");
    exit(1);
  }
  for (int i = 0; i < size; i++) {
    pixels[i] = malloc(sizeof(struct ppm_pixel) * size);
    if (pixels[i] == NULL) {
      printf("Error: failed malloc.  Exiting...\n");
      exit(1);
    }
  }

  // allocate memory for array of counts
  counts = malloc(sizeof(struct pix_counts *) * size);
  if (counts == NULL) {
    printf("Error: failed malloc.  Exiting...\n");
    exit(1);
  }
  for (int i = 0; i < size; i++) {
    counts[i] = malloc(sizeof(struct pix_counts) * size);
    if (counts[i] == NULL) {
      printf("Error: failed malloc.  Exiting...\n");
      exit(1);
    }
  }
  // set all struct variables to zero
  initializeCounts(counts, size);

  // allocate memory for palette and fill with colors
  palette = malloc(sizeof(struct ppm_pixel) * 256);
  if (palette == NULL) {
    printf("Error: failed malloc.  Exiting...\n");
    exit(1);
  }
  fillPalette(palette, "palette30.txt");

  gettimeofday(&tstart, NULL);

  // pick random point as seed of orbit, with x,y in range [-1,1]
  struct point seed;
  seed.x = 2 * ((float) rand() / RAND_MAX) - 1;
  seed.y = 2 * ((float) rand() / RAND_MAX) - 1;
  // random initial color with r,b,g in [0, 1]
  float initR = (float) rand() / RAND_MAX;
  float initG = (float) rand() / RAND_MAX;
  float initB = (float) rand() / RAND_MAX;
  iterate(size, xmin, xmax, ymin, ymax, seed, initR, initG, initB, iterations,
      counts, palette, &maxCount);
  computeColor(pixels, counts, palette, &maxCount, size);







  gettimeofday(&tend, NULL);
  timer = tend.tv_sec - tstart.tv_sec + (tend.tv_usec - tstart.tv_usec)/1.e6;
  printf("Computed fractal flame (%dx%d) in %g seconds\n", size, size, timer);

  // write to file
  new_file[0] = '\0';
  sprintf(new_file, "flame_S%d_%lu.ppm", size, time(0));
  printf("Writing file %s\n", new_file);
  write_ppm(new_file, pixels, size, size);

  // free allocated array memory
  for (int i = 0; i < size; i++) {
    free(pixels[i]);
    pixels[i] = NULL;
    free(counts[i]);
    counts[i] = NULL;
  }
  free(pixels);
  pixels = NULL;
  free(counts);
  counts = NULL;
  free(palette);
  palette = NULL;

  return 0;
}
