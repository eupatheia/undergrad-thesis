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
#include "quicksort.h"

// array of function pointers
transform_ptr transforms[17] = {linear, sinusoidal, spherical, swirl,
    horseshoe, polar, handkerchief, heart, disk, spiral, hyperbolic, diamond,
    ex, julia, fisheye, exponential, eyefish};

char * funcNames[17] = {"linear", "sinusoidal", "spherical", "swirl",
    "horseshoe", "polar", "handkerchief", "heart", "disk", "spiral",
    "hyperbolic", "diamond", "ex", "julia", "fisheye", "exponential",
    "eyefish"};

// save second generated Gaussian for later
float noise = 0;  // nothing yet
int savedNoise = 0;  // if 1, use noise instead of regenerating

// print iterated system info to console
void printSystemStats(struct systemInfo * info) {
  printf("F0: weight = %.3f\n(%.6f, %.6f, %.6f, %.6f, %.6f, %.6f)\n",
      info->weights[0], info->affineParams[0], info->affineParams[1],
      info->affineParams[2], info->affineParams[3], info->affineParams[4],
      info->affineParams[5]);
  printf("F1: weight = %.3f\n(%.6f, %.6f, %.6f, %.6f, %.6f, %.6f)\n",
      info->weights[1], info->affineParams[6], info->affineParams[7],
      info->affineParams[8], info->affineParams[9], info->affineParams[10],
      info->affineParams[11]);
  if (info->numFunctions > 2) {
    printf("F2: weight = %.3f\n(%.6f, %.6f, %.6f, %.6f, %.6f, %.6f)\n",
      info->weights[2], info->affineParams[12], info->affineParams[13],
      info->affineParams[14], info->affineParams[15], info->affineParams[16],
      info->affineParams[17]);
  }
  if (info->numFunctions > 3) {
    printf("F2: weight = %.3f\n(%.6f, %.6f, %.6f, %.6f, %.6f, %.6f)\n",
      info->weights[2], info->affineParams[18], info->affineParams[19],
      info->affineParams[20], info->affineParams[21], info->affineParams[22],
      info->affineParams[23]);
  }
}

// generate new random function system with
//   2 to maxFunctions functions, where maxFunctions = {2,3,4}, ideally
//   1 to maxSymmetry symmetry, where maxSymmetry = {1,2,3}
//   6 affine parameters per function
struct systemInfo getSystem(int maxSymmetry, int maxFunctions) {
  struct systemInfo info;
  info.numFunctions = (rand() % (maxFunctions - 1)) + 2;
  info.symmetry = (rand() % maxSymmetry) + 1;
  info.functions = malloc(sizeof(transform_ptr) * info.numFunctions);
  info.weights = malloc(sizeof(float) * info.numFunctions);
  int numAffine = info.numFunctions * 6;
  info.affineParams = malloc(sizeof(float) * numAffine);

  printf("Num functions: %d\n", info.numFunctions);
  printf("Symmetry: %d\n", info.symmetry);

  for (int i = 0; i < info.numFunctions; i++) {
    int choice = rand() % 17;
    info.functions[i] = transforms[choice];
    printf("F%d: %s\n", i, funcNames[choice]);
  }
  float prob;  // function weights must sum to this (excluding rotations)
  if (info.symmetry == 1) {
    prob = 1;
  } else if (info.symmetry == 2) {
    prob = 0.5;
  } else {
    info.symmetry = 3;  // clamp b/c not handling higher symmetry
    prob = 1.0 / 3.0;
  }
  for (int i = 0; i < info.numFunctions - 1; i++) {
    info.weights[i] = randomParam(0, prob);
    prob -= info.weights[i];
  }
  info.weights[info.numFunctions - 1] = prob;
  for (int i = 0; i < numAffine; i++) {
    info.affineParams[i] = randomParam(-1, 1);
  }
  printSystemStats(&info);
  return info;
}

// Iterate through the function system and accumulate counts,
// return quality value
void iterate(struct systemInfo * info, int size, float xmin, float xmax,
    float ymin, float ymax, struct point p, float c, int iterations,
    struct pix_counts ** counts, struct ppm_pixel * palette, int * maxCount) {

  int yrow, xcol;
  void (*systemToIterate) (struct point *, float *, struct systemInfo *);
  if (info->symmetry == 1) {
    systemToIterate = system1Sym;
  } else if (info->symmetry == 2) {
    systemToIterate = system2Sym;
  } else {
    systemToIterate = system3Sym;
  }
  for (int i = 0; i < iterations; i++) {
    // pick from a system of functions and calculate new point and color
    systemToIterate(&p, &c, info);
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
void renderSupersample(int start_row, int end_row, int start_col, int end_col,
    struct ppm_pixel ** pixels, struct pix_counts ** counts,
    struct ppm_pixel * palette, int internalSize, int outputSize,
    int * maxCount, float * gaussian3, float * gaussian5, float * gaussian7) {
  for (int i = start_row; i < end_row; i++) {
    for (int j = start_col; j < end_col; j++) {
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

// helper function to calculate quality of a specimen in terms of
// proportion of black pixels (never reached)
float calculateQuality(struct pix_counts ** counts, int internalSize) {
  int colorful = 0;  // num of non-black indices 
  for (int i = 0; i < internalSize; i++) {
    for (int j = 0; j < internalSize; j++) {
      if (counts[i][j].alpha > 0) {
        colorful++;
      }
    }
  }
  return (float) colorful / (internalSize * internalSize);
}

// helper function to sample from a Guassian distribution N(mean, stddev),
// given pseudo-random numbers from a uniform distribution,
// using the Box-Mueller Transform:
// see https://en.wikipedia.org/wiki/Box–Muller_transform
float randGaussian(float mean, float stddev) {
  if (savedNoise) {
    savedNoise = 0;  // reset
    return noise;
  }
  float rand1 = randomParam(0, 1);
  float rand2 = randomParam(0, 1);
  float radius = sqrt(-2 * log(rand1));
  float z = stddev * radius * cos(2 * M_PI * rand2) + mean;
  noise = stddev * radius * sin(2 * M_PI * rand2) + mean;  // saved for later
  savedNoise = 1;  // next time use saved noise
  return z;
}

// helper function makes a copy of parent into child, then mutates weights and
// affine params by some random number sampled from a Gaussian distribution,
// assuming functions, weights, and affineParams have been freed
void mutate(struct systemInfo * child, struct systemInfo * parent) {
  child->numFunctions = parent->numFunctions;
  child->symmetry = parent->symmetry;
  child->functions = malloc(sizeof(transform_ptr) * child->numFunctions);
  child->weights = malloc(sizeof(float) * child->numFunctions);
  int numAffine = child->numFunctions * 6;
  child->affineParams = malloc(sizeof(float) * numAffine);

  // copy over functions
  for (int i = 0; i < child->numFunctions; i++) {
    child->functions[i] = parent->functions[i];
  }

  // change params by some small random number
  for (int i = 0; i < child->numFunctions; i++) {
    child->weights[i] = parent->weights[i] + randGaussian(0, 0.1);
  }
  for (int i = 0; i < numAffine; i++) {
    child->affineParams[i] = parent->affineParams[i] + randGaussian(0, 0.1);
  }
}

// helper function to free system info
void freeSystemInfo(struct systemInfo * info) {
  // free system info
  free(info->functions);
  info->functions = NULL;
  free(info->weights);
  info->weights = NULL;
  free(info->affineParams);
  info->affineParams = NULL;
}

// helper function to calculate the max of two ints
int max(int a, int b) {
  if (a > b) {
    return a;
  } else {
    return b;
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

// helper function to divide plane into sections and assign coordinates,
// where each section is like a "stripe" across the plane,
// and i is a section number from 0 to numSections - 1
void getCoordinates(int i, int numSections, int size, int * start_col,
    int * end_col, int * start_row, int * end_row) {
  // for now, only allow even division of the plane
  if (size % numSections != 0) {
    printf("error: %d does not divide evenly by %d\n", size, numSections);
    exit(1);
  }
  *start_col = 0;
  *end_col = size;
  *start_row = (size / numSections) * i;
  *end_row = (size / numSections) * (i + 1);
}

struct thread_data {
  int id;
  int start_row;
  int end_row;
  int start_col;
  int end_col;
  struct ppm_pixel ** pixels;
  struct pix_counts ** counts;
  struct ppm_pixel * palette;
  int internalSize;
  int outputSize;
  int * maxCount;
  float * gaussian3;
  float * gaussian5;
  float * gaussian7;
};

void * thread_function(void * args) {
  struct thread_data * data = (struct thread_data *) args;
  // printf("Thread %d) sub-image block: cols (%d, %d) to rows (%d, %d)\n",
  //     data->id, data->start_col, data->end_col, data->start_row, data->end_row);

  renderSupersample(data->start_col, data->end_col, data->start_row,
    data->end_row, data->pixels, data->counts, data->palette,
    data->internalSize, data->outputSize, data->maxCount, data->gaussian3,
    data->gaussian5, data->gaussian7);

  // printf("Thread %d) finished\n", data->id);
  return (void *) NULL;
}

// helper function to write a particular specimen to a file
void threadOutputSpecimen(int internalSize, int outputSize, int iterations,
    int numProcesses, struct ppm_pixel ** pixels, struct ppm_pixel * palette,
    int paletteNum, pthread_t * threads, struct thread_data * data,
    float * gaussian3, float * gaussian5, float * gaussian7,
    struct pix_counts ** counts, struct timeval tstart, int * maxCount,
    int doSearch, int generation) {

  int start_col, end_col, start_row, end_row;
  char new_file[128];
  struct timeval tend;
  double timer;

  for (int i = 0; i < numProcesses; i++) {
    getCoordinates(i, numProcesses, outputSize, &start_col, &end_col,
        &start_row, &end_row);
    data[i].id = i;
    data[i].start_row = start_row;
    data[i].end_row = end_row;
    data[i].start_col = start_col;
    data[i].end_col = end_col;
    data[i].pixels = pixels;
    data[i].counts = counts;
    data[i].palette = palette;
    data[i].internalSize = internalSize;
    data[i].outputSize = outputSize;
    data[i].maxCount = maxCount;
    data[i].gaussian3 = gaussian3;
    data[i].gaussian5 = gaussian5;
    data[i].gaussian7 = gaussian7;
    // create threads
    pthread_create(&threads[i], NULL, thread_function, (void *) &data[i]);
  }

  // join threads
  for (int i = 0; i < numProcesses; i++) {
    pthread_join(threads[i], NULL);
  }

  gettimeofday(&tend, NULL);
  timer = tend.tv_sec - tstart.tv_sec + (tend.tv_usec - tstart.tv_usec)/1.e6;
  printf("Computed fractal flame (%dx%d) in %g seconds\n",
      outputSize, outputSize, timer);

  // write to file
  new_file[0] = '\0';
  if (doSearch) {
    sprintf(new_file, "gen%d_flame_S%d_N%d_C%d_%lu.ppm", generation, outputSize,
      iterations, paletteNum, time(0));
  } else {
    sprintf(new_file, "flame_S%d_N%d_C%d_%lu.ppm", outputSize, iterations,
      paletteNum, time(0));
  }
  printf("Writing file %s\n\n", new_file);
  write_ppm(new_file, pixels, outputSize, outputSize);
}

// generates a specific flame by iteration and multithread coloring;
// if doSearch = 1 = true, do mu + lambda evolutionary search,
// else randomly generate one flame
void generate(int internalSize, int outputSize, float xmin, float xmax,
    float ymin, float ymax, int iterations, int numProcesses,
    struct ppm_pixel ** pixels, struct ppm_pixel * palette, int paletteNum,
    pthread_t * threads, struct thread_data * data, float * gaussian3,
    float * gaussian5, float * gaussian7, int doSearch) {

  struct timeval tstart;
  int mu = 1;  // num retained per generation
  int lambda = 1;  // num reproduced per generation
  int populationSize;
  // array of 2D count arrays represents the flames in one generation
  struct specimen * population;
  int generations = 0;

  if (doSearch) {
    populationSize = mu + lambda;
  } else {
    populationSize = 1;
  }

  // allocate memory for population
  population = malloc(sizeof(struct specimen) * (populationSize));
  if (population == NULL) {
    printf("Error: failed malloc.  Exiting...\n");
    exit(1);
  }
  for (int i = 0; i < populationSize; i++) {
    population[i].counts = malloc(sizeof(struct pix_counts *) * internalSize);
    if (population[i].counts == NULL) {
      printf("Error: failed malloc.  Exiting...\n");
      exit(1);
    }
    for (int j = 0; j < internalSize; j++) {
      population[i].counts[j] = malloc(sizeof(struct pix_counts) * internalSize);
      if (population[i].counts[j] == NULL) {
        printf("Error: failed malloc.  Exiting...\n");
        exit(1);
      }
    }
  }

  // pick random point as seed of orbit, with x,y in range [-1,1]
  struct point seed;
  seed.x = 2 * ((float) rand() / RAND_MAX) - 1;
  seed.y = 2 * ((float) rand() / RAND_MAX) - 1;
  // random initial color in [0, 1]
  float initColor = (float) rand() / RAND_MAX;

  gettimeofday(&tstart, NULL);

  if (doSearch) {
    // random initial specimens
    for (int i = 0; i < populationSize; i++) {
      // set all counts to zero
      initializeCounts(population[i].counts, internalSize);
      population[i].maxCount = 0;
      population[i].info = getSystem(1, 4);
      iterate(&(population[i].info), internalSize, xmin, xmax, ymin, ymax, seed,
          initColor, iterations, population[i].counts, palette,
          &(population[i].maxCount));
      population[i].quality = calculateQuality(population[i].counts, internalSize);

      // write to file
      printf("Quality %.3f\n", population[i].quality);
      threadOutputSpecimen(internalSize, outputSize, iterations, numProcesses,
          pixels, palette, paletteNum, threads, data, gaussian3, gaussian5,
          gaussian7, population[i].counts, tstart, &population[i].maxCount,
          doSearch, generations);
    }
    generations++;
    hybridSort(population, 0, populationSize - 1);
    // do evolutionary search find a high quality specimen (quality >= 0.5)
    while (population[populationSize - 1].quality < 0.3) {
      // start another generation
      // replace the lambda lowest quality specimens
      for (int i = 0; i < lambda; i++) {
        freeSystemInfo(&population[i].info);
        // set all counts to zero again
        initializeCounts(population[i].counts, internalSize);
        population[i].maxCount = 0;
        // make mutated copies of higher-quality specimens
        mutate(&(population[i].info), &(population[max(i + lambda, populationSize - 1)].info));
        iterate(&(population[i].info), internalSize, xmin, xmax, ymin, ymax,
            seed, initColor, iterations, population[i].counts, palette,
            &(population[i].maxCount));
        population[i].quality = calculateQuality(population[i].counts,
            internalSize);

        printf("Quality %.3f\n", population[i].quality);
        // write to file using multiple threads
        threadOutputSpecimen(internalSize, outputSize, iterations, numProcesses,
            pixels, palette, paletteNum, threads, data, gaussian3, gaussian5,
            gaussian7, population[i].counts, tstart, &population[i].maxCount,
            doSearch, generations);
      }
      generations++;
      hybridSort(population, 0, populationSize - 1);
    }
  } else {
    // set all counts to zero
    initializeCounts(population[0].counts, internalSize);
    population[0].maxCount = 0;
    population[0].info = getSystem(3, 4);
    // output one random flame
    iterate(&(population[0].info), internalSize, xmin, xmax, ymin, ymax, seed,
        initColor, iterations, population[0].counts, palette,
        &(population[0].maxCount));
    population[0].quality = calculateQuality(population[0].counts, internalSize);
  }

  // print out the highest quality specimen (or the one random specimen)
  struct specimen bestSpecimen = population[populationSize - 1];
  printf("%d generations to get final quality %.3f\n", generations,
      bestSpecimen.quality);

  // write to file using multiple threads
  threadOutputSpecimen(internalSize, outputSize, iterations, numProcesses,
    pixels, palette, paletteNum, threads, data, gaussian3, gaussian5, gaussian7,
    bestSpecimen.counts, tstart, &bestSpecimen.maxCount, doSearch, generations);

  // free arrays
  for (int i = 0; i < populationSize; i++) {
    freeSystemInfo(&population[i].info);
  }
  for (int i = 0; i < populationSize; i++) {
    for (int j = 0; j < internalSize; j++) {
      free(population[i].counts[j]);
      population[i].counts[j] = NULL;
    }
    free(population[i].counts);
    population[i].counts = NULL;
  }
  free(population);
  population = NULL;
}

//===========================================================//
//===========================================================//
//===========================================================//

int main(int argc, char* argv[]) {
  srand(time(0));  // give random seed to generator
  int palettes[12] = {11, 13, 20, 26, 30, 31, 34, 35, 41, 57, 82, 87};

  int outputSize = 900;
  int internalSize = outputSize * 3;
  // use bi-unit square
  float xmin = -1.5;
  float xmax = 1.5;
  float ymin = -1.5;
  float ymax = 1.5;
  int iterations = 40000000;
  int paletteNum = palettes[rand() % 12];
  int numProcesses = 4;
  struct ppm_pixel ** pixels = NULL;
  struct ppm_pixel * palette = NULL;
  char new_file[100];
  pthread_t * threads;
  struct thread_data * data;
  // global normalized gaussian blur matrices
  float * gaussian3;
  float * gaussian5;
  float * gaussian7;
  int doSearch = 0;  // default no search, i.e. generate one randomly

  int opt;
  while ((opt = getopt(argc, argv, ":s:n:c:p:e")) != -1) {
    switch (opt) {
      case 's': outputSize = atoi(optarg);
                internalSize = outputSize * 3; break;
      case 'n': iterations = atoi(optarg); break;
      case 'c': paletteNum = atoi(optarg); break;
      case 'p': numProcesses = atoi(optarg); break;
      case 'e': doSearch = 1; break;  // no colon in opt string b/c no arg
      case '?': printf("usage: %s -s <outputSize> -n <iterations>"
          " -c <paletteNumber> -p <numProcesses> -e [no args, do evolutionary"
          "search instead of random generation]\n", argv[0]); break;
    }
  }
  printf("Generating flame with size %dx%d\n", outputSize, outputSize);
  printf("  Iterations = %d\n", iterations);
  printf("  Color Palette = %d\n", paletteNum);
  printf("  Num processes = %d\n", numProcesses);

  // allocate memory for thread identifiers
  threads = malloc(sizeof(pthread_t) * numProcesses);
  if (threads == NULL) {
    printf("Error: failed malloc.  Exiting...\n");
    exit(1);
  }
  // allocate memory for thread function data structs
  data = malloc(sizeof(struct thread_data) * numProcesses);
  if (data == NULL) {
    printf("Error: failed malloc.  Exiting...\n");
    exit(1);
  }

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

  // allocate memory for palette and fill with colors
  palette = malloc(sizeof(struct ppm_pixel) * 256);
  if (palette == NULL) {
    printf("Error: failed malloc.  Exiting...\n");
    exit(1);
  }
  new_file[0] = '\0';
  sprintf(new_file, "palettes/palette%d.txt", paletteNum);
  fillPalette(palette, new_file);

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

  generate(internalSize, outputSize, xmin, xmax, ymin, ymax, iterations,
      numProcesses, pixels, palette, paletteNum, threads, data,
      gaussian3, gaussian5, gaussian7, doSearch);

  // free allocated array memory
  for (int i = 0; i < outputSize; i++) {
    free(pixels[i]);
    pixels[i] = NULL;
  }
  free(pixels);
  pixels = NULL;
  free(palette);
  palette = NULL;
  free(gaussian3);
  gaussian3 = NULL;
  free(gaussian5);
  gaussian5 = NULL;
  free(gaussian7);
  gaussian7 = NULL;
  free(threads);
  threads = NULL;
  free(data);
  data = NULL;

  return 0;
}
