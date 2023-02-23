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
      counts[yrow][xcol].color = c;  // set blended color
      counts[yrow][xcol].alpha++;  // increment counter
      // update max values
      if (counts[yrow][xcol].alpha > *maxCount) {
        *maxCount = counts[yrow][xcol].alpha;
      }
    }
  }
}

// no supersampling, white if reached, black otherwise
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

// no supersampling, grayscale color based on log frequency of hits
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

// color by supersampling and reduce image to outputSize
void renderSupersample(struct ppm_pixel ** pixels, struct pix_counts ** counts,
    struct ppm_pixel * palette, int internalSize, int outputSize,
    int * maxCount, float * gaussian3, float * gaussian5, float * gaussian7) {
  for (int i = 0; i < outputSize; i++) {
    for (int j = 0; j < outputSize; j++) {
      int oldi = i * 3 + 1;  // position in internal array
      int oldj = j * 3 + 1;  // position in internal array
      // diameter (num pixels from each side to center, not including center)
      // inversely proportional to hit count, so range 1-3 (i.e. max 7x7)
      int diameter;
      if (counts[oldi][oldj].alpha == 0) {
        diameter = 1;
      } else if (counts[oldi][oldj].alpha == 1) {
        diameter = 3;
      } else {
        diameter = 3 * (1.0 / counts[oldi][oldj].alpha) + 1;
      }

      float * kernel;
      if (diameter == 1) {
        kernel = gaussian3;
      } else if (diameter == 2) {
        kernel = gaussian5;
      } else if (diameter == 3) {
        kernel = gaussian7;
      } else {
        // logical error, stop computing and return
        printf("ERROR: kernel diameter is %d (not 1, 2, or 3)\n", diameter);
        return;
      }

      float countAvg = 0.0;
      float rAvg = 0.0;
      float gAvg = 0.0;
      float bAvg = 0.0;
      int width = diameter * 2 + 1;
      int samples = 0;
      int startRow = fmax(oldi - diameter, 0);
      int endRow = fmin(oldi + diameter, internalSize - 1);
      int startCol = fmax(oldj - diameter, 0);
      int endCol = fmin(oldj + diameter, internalSize - 1);
      for (int m = startRow; m <= endRow; m++) {
        for (int n = startCol; n <= endCol; n++) {
          int ki = m - (oldi) + diameter;
          int kj = n - (oldj) + diameter;
          float kij = kernel[ki * width + kj];
          countAvg += counts[m][n].alpha * kij;
          struct ppm_pixel origColor = palette[(int) (counts[m][n].color * 255)];
          rAvg += origColor.red * kij;
          gAvg += origColor.green * kij;
          bAvg += origColor.blue * kij;
          samples++;
        }
      }
      float factor;
      if (countAvg < 1) {
        factor = 0.000001;
      } else {
        factor = log(countAvg) / log(*maxCount);
      }
      int brightness = 5;
      float gamma = 1;
      pixels[i][j].red = fmin(pow((rAvg / 255.0) * factor, 1.0 / gamma)
          * brightness * 255, 255);
      pixels[i][j].green = fmin(pow((gAvg / 255.0) * factor, 1.0 / gamma)
          * brightness * 255, 255);
      pixels[i][j].blue = fmin(pow((bAvg / 255.0) * factor, 1.0 / gamma)
          * brightness * 255, 255);
    }
  }
}

// no supersampling, assuming pixels and counts are same size
void renderNoSupersample(struct ppm_pixel ** pixels,
    struct pix_counts ** counts, struct ppm_pixel * palette, int internalSize,
    int outputSize, int * maxCount) {
  // check that sizes match
  if (internalSize != outputSize) {
    printf("ERROR: cannot render if internalSize != outputSize\n");
    return;
  }
  for (int i = 0; i < outputSize; i++) {
    for (int j = 0; j < outputSize; j++) {
      if (counts[i][j].alpha > 0) {
        float factor = log(counts[i][j].alpha) / log(*maxCount);
        int index = counts[i][j].color * 255;
        int brightness = 1;
        float gamma = 1;
        pixels[i][j].red = fmin(pow((palette[index].red / 255.0) * factor, 
            1.0 / gamma) * brightness * 255, 255);
        pixels[i][j].green =  fmin(pow((palette[index].green / 255.0) * factor,
            1.0 / gamma) * brightness * 255, 255);
        pixels[i][j].blue =  fmin(pow((palette[index].blue / 255.0) * factor,
            1.0 / gamma) * brightness * 255, 255);
      } else {
        // never reached, color black
        pixels[i][j].red = 0;
        pixels[i][j].green = 0;
        pixels[i][j].blue = 0;
      }
    }
  }
}


// helper function to calculate normalized gaussian blur matrix
void calculateBlurMatrix(float * matrix, int diameter) {
  float sum = 0.0;
  for (int m = -diameter; m <= diameter; m++) {
    for (int n = -diameter; n <= diameter; n++) {
      int index = matrix[((diameter + m) * (diameter * 2 + 1)) + (diameter + n)];
      float val = (1.0 / (2 * M_PI)) * exp(-(pow(m, 2) + pow(n, 2)) / 2.0);
      sum += val;
      matrix[((diameter + m) * (diameter * 2 + 1)) + (diameter + n)] = val;
    }
  }
  // normalize
  for (int m = -diameter; m <= diameter; m++) {
    for (int n = -diameter; n <= diameter; n++) {
      matrix[((diameter + m) * (diameter * 2 + 1)) + (diameter + n)] /= sum;
    }
  }
}

// helper function to set all struct variables to zero
// for ease of incrementation
void initializeCounts(struct pix_counts ** counts, int size) {
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      counts[i][j].color = 0.0;
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
  int outputSize = 900;
  int internalSize = outputSize * 3;
  // use bi-unit square
  float xmin = -1.5;
  float xmax = 1.5;
  float ymin = -1.5;
  float ymax = 1.5;
  int iterations = 40000000;
  struct ppm_pixel ** pixels = NULL;
  struct pix_counts ** counts = NULL;
  struct ppm_pixel * palette = NULL;
  double timer;
  struct timeval tstart, tend;
  char new_file[100];
  float x, y;
  int k;  // to randomly choose a function
  int maxCount = 0;
  // global normalized gaussian blur matrices
  float * gaussian3;
  float * gaussian5;
  float * gaussian7;

  srand(time(0));  // give random seed to generator

  int opt;
  while ((opt = getopt(argc, argv, ":s:n:")) != -1) {
    switch (opt) {
      case 's': outputSize = atoi(optarg);
                internalSize = outputSize * 3; break;
      case 'n': iterations = atoi(optarg); break;
      case '?': printf("usage: %s -s <outputSize> -n <iterations>\n", argv[0]); break;
    }
  }
  printf("Generating flame with size %dx%d\n", outputSize, outputSize);
  printf("    and %d iterations\n", iterations);

  // allocate memory for output pixels
  pixels = malloc(sizeof(struct ppm_pixel *) * outputSize);
  if (pixels == NULL) {
    printf("Error: failed malloc.  Exiting...\n");
    exit(1);
  }
  for (int i = 0; i < outputSize; i++) {
    pixels[i] = malloc(sizeof(struct ppm_pixel) * outputSize);
    if (pixels[i] == NULL) {
      printf("Error: failed malloc.  Exiting...\n");
      exit(1);
    }
  }

  // allocate memory for internal array of counts
  counts = malloc(sizeof(struct pix_counts *) * internalSize);
  if (counts == NULL) {
    printf("Error: failed malloc.  Exiting...\n");
    exit(1);
  }
  for (int i = 0; i < internalSize; i++) {
    counts[i] = malloc(sizeof(struct pix_counts) * internalSize);
    if (counts[i] == NULL) {
      printf("Error: failed malloc.  Exiting...\n");
      exit(1);
    }
  }
  // set all struct variables to zero
  initializeCounts(counts, internalSize);

  // allocate memory for palette and fill with colors
  palette = malloc(sizeof(struct ppm_pixel) * 256);
  if (palette == NULL) {
    printf("Error: failed malloc.  Exiting...\n");
    exit(1);
  }
  fillPalette(palette, "palette30.txt");

  // allocate memory for kernels
  gaussian3 = malloc(sizeof(float) * 3 * 3);
  if (gaussian3 == NULL) {
    printf("Error: failed malloc.  Exiting...\n");
    exit(1);
  }
  gaussian5 = malloc(sizeof(float) * 5 * 5);
  if (gaussian5 == NULL) {
    printf("Error: failed malloc.  Exiting...\n");
    exit(1);
  }
  gaussian7 = malloc(sizeof(float) * 7 * 7);
  if (gaussian7 == NULL) {
    printf("Error: failed malloc.  Exiting...\n");
    exit(1);
  }

  // calculate normalized Gaussian blur matrix
  calculateBlurMatrix(gaussian3, 1);
  calculateBlurMatrix(gaussian5, 2);
  calculateBlurMatrix(gaussian7, 3);

  // pick random point as seed of orbit, with x,y in range [-1,1]
  struct point seed;
  seed.x = 2 * ((float) rand() / RAND_MAX) - 1;
  seed.y = 2 * ((float) rand() / RAND_MAX) - 1;
  // random initial color in [0, 1]
  float initColor = (float) rand() / RAND_MAX;

  gettimeofday(&tstart, NULL);

  iterate(internalSize, xmin, xmax, ymin, ymax, seed, initColor, iterations,
      counts, palette, &maxCount);
  renderSupersample(pixels, counts, palette, internalSize, outputSize,
      &maxCount,gaussian3, gaussian5, gaussian7);
  // renderNoSupersample(pixels, counts, palette, internalSize, outputSize,
  //     &maxCount);



  
  







  gettimeofday(&tend, NULL);
  timer = tend.tv_sec - tstart.tv_sec + (tend.tv_usec - tstart.tv_usec)/1.e6;
  printf("Computed fractal flame (%dx%d) in %g seconds\n",
      outputSize, outputSize, timer);

  // write to file
  new_file[0] = '\0';
  sprintf(new_file, "flame_S%d_N%d_%lu.ppm", outputSize, iterations, time(0));
  printf("Writing file %s\n", new_file);
  write_ppm(new_file, pixels, outputSize, outputSize);

  // free allocated array memory
  for (int i = 0; i < outputSize; i++) {
    free(pixels[i]);
    pixels[i] = NULL;
  }
  for (int i = 0; i < internalSize; i++) {
    free(counts[i]);
    counts[i] = NULL;
  }
  free(pixels);
  pixels = NULL;
  free(counts);
  counts = NULL;
  free(palette);
  palette = NULL;
  free(gaussian3);
  gaussian3 = NULL;
  free(gaussian5);
  gaussian5 = NULL;
  free(gaussian7);
  gaussian7 = NULL;

  return 0;
}
