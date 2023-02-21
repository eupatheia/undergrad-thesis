#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <assert.h>
#include <time.h>
#include <pthread.h>
#include <sys/time.h>
#include <math.h>
#include "transforms.h"
#include "read_write.h"

// Iterate through the function system and accumulate counts
void iterate(int size, float xmin, float xmax, float ymin, float ymax, 
    struct point p, float c, int iterations, struct pix_counts ** counts,
    struct ppm_pixel * palette, int * maxCount) {

  // associated color for F_0
  float c0 = 1.0;
  // associated color for F_1
  float c1 = 0.0;

  int k, yrow, xcol;
  for (int i = 0; i < iterations; i++) {
    k = rand() % 2;  // both have probability 0.5, here
    if (k == 0) {
      p = affine(p, 0.562482, -0.539599, -0.42992, 0.397861, 0.501088, -0.112404);
      c = (c + c0) / 2;
    } else {  // k == 1
      p = affine(p, 0.830039, 0.16248, 0.91022, -0.496174, 0.75046, 0.288389);
      c = (c + c1) / 2;
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
      counts[yrow][xcol].sum = c;
      counts[yrow][xcol].alpha++;
      // update max values
      if (counts[yrow][xcol].alpha > *maxCount) {
        *maxCount = counts[yrow][xcol].alpha;
      }
    }
  }
}

void computeColorBinary(struct ppm_pixel ** pixels, struct pix_counts ** counts,
    int * maxCount, int size) {
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      if (counts[i][j].alpha > 0) {
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
      int intensity = (log(counts[i][j].alpha) / log(*maxCount)) * 255;
      pixels[i][j].red = intensity;
      pixels[i][j].green = intensity;
      pixels[i][j].blue = intensity;
    }
  }
}

void computeColor(struct ppm_pixel ** pixels, struct pix_counts ** counts,
    struct ppm_pixel * palette, int size) {
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      if (counts[i][j].alpha > 0) {
        float factor = log(counts[i][j].alpha) / counts[i][j].alpha;
        int index = counts[i][j].sum * 255;
        pixels[i][j] = palette[index];
        // pixels[i][j].red = round(pow((palette[index].red / 255.0f), factor) * 255.0f);
        // pixels[i][j].green = round(pow((palette[index].green / 255.0f), factor) * 255.0f);
        // pixels[i][j].blue = round(pow((palette[index].blue / 255.0f), factor) * 255.0f);
        printf("%d\n", index);
      } else {
        // never visited, color black
        pixels[i][j].red = 0;
        pixels[i][j].green = 0;
        pixels[i][j].blue = 0;
      }
    }
  }
}

// helper function to set all struct variables to zero
// for ease of incrementation
void initializeCounts(struct pix_counts ** counts, int size) {
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      counts[i][j].sum = 0.0;
      counts[i][j].alpha = 0;
    }
  }
}

// helper function to fill 16x16 ppm file with palette colors
void outputPalette(struct ppm_pixel ** pixels, struct ppm_pixel * palette,
    int size) {
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      pixels[i][j] = palette[i * 16 + j];
    }
  }
}

//===========================================================//
//===========================================================//
//===========================================================//

int main(int argc, char* argv[]) {
  int size = 900;
  // use bi-unit square
  float xmin = -1.5;
  float xmax = 1.5;
  float ymin = -1.5;
  float ymax = 1.5;
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

  srand(time(0));  // give random seed to generator

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
  // random initial color in [0, 1]
  float initColor = (float) rand() / RAND_MAX;
  iterate(size, xmin, xmax, ymin, ymax, seed, initColor, iterations,
      counts, palette, &maxCount);
  computeColor(pixels, counts, palette, size);



  
  







  gettimeofday(&tend, NULL);
  timer = tend.tv_sec - tstart.tv_sec + (tend.tv_usec - tstart.tv_usec)/1.e6;
  printf("Computed fractal flame (%dx%d) in %g seconds\n", size, size, timer);

  // write to file
  new_file[0] = '\0';
  sprintf(new_file, "flame_S%d_N%d_%lu.ppm", size, iterations, time(0));
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
