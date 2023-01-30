#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <assert.h>
#include <time.h>
#include <pthread.h>
#include <sys/time.h>
#include <math.h>
#include "read_ppm.h"

pthread_barrier_t barrier;
pthread_mutex_t mutex;

// Step 1: Determine mandelbrot set membership
void computeMembership(int size, int start_row, int end_row, int start_col,
    int end_col, float xmin, float xmax, float ymin, float ymax,
    int maxIterations, int ** in_set) {

  float xfrac, yfrac, x0, y0, x, y, xtmp;
  int iter;

  for (int i = start_row; i < end_row; i++) {
    for (int j = start_col; j < end_col; j++) {
      xfrac = (float) j / size;
      yfrac = (float) i / size;
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
        in_set[i][j] = 0;
      } else {
        // did not escape, in set
        in_set[i][j] = 1;
      }
    }
  }
}

// Step 2: Compute visited counts
void computeCounts(int size, int start_row, int end_row, int start_col,
    int end_col, float xmin, float xmax, float ymin, float ymax,
    int * max_count, int ** in_set, int ** counts) {

  float xfrac, yfrac, x0, y0, x, y, xtmp;
  int yrow, xcol;

  for (int i = start_row; i < end_row; i++) {
    for (int j = start_col; j < end_col; j++) {
      if (in_set[i][j] == 1) {
        continue;
      } else {
        xfrac = (float) j / size;
        yfrac = (float) i / size;
        x0 = xmin + xfrac * (xmax - xmin);
        y0 = ymin + yfrac * (ymax - ymin);

        x = 0;
        y = 0;
        while (x*x + y*y < 2*2) {
          xtmp = x*x - y*y + x0;
          y = 2*x*y + y0;
          x = xtmp;

          yrow = round(size * (y - ymin) / (ymax - ymin));
          xcol = round(size * (x - xmin) / (xmax - xmin));
          if (yrow < 0 || yrow >= size) {
            continue; // out of range
          }
          if (xcol < 0 || xcol >= size) {
            continue; // out of range
          }

	        pthread_mutex_lock(&mutex);
          counts[yrow][xcol]++;
          // update max count
          if (counts[yrow][xcol] > *max_count) {
            *max_count = counts[yrow][xcol];
          }
          pthread_mutex_unlock(&mutex);
	      }
      }
    }
  }
}

// Step 3: Compute colors
void computeColors(int start_row, int end_row, int start_col, int end_col,
    int * max_count, int ** counts, struct ppm_pixel ** pixels, int channel) {

  float gamma = 0.681;
  float factor = 1.0 / gamma;

  for (int i = start_row; i < end_row; i++) {
    for (int j = start_col; j < end_col; j++) {
      float value = 0;
      if (counts[i][j] > 0) {
	      value = log(counts[i][j]) / log(*max_count);
        value = pow(value, factor);
      }
      if (channel == 1) {
        pixels[i][j].red = value * 255;
      } else if (channel == 2) {
        pixels[i][j].green = value * 255;
      } else { // channel == 3
        pixels[i][j].blue = value * 255;
      }
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
    printf("error: %d does not divide evenly by %d", size, numSections);
    exit(1);
  }
  *start_col = 0;
  *end_col = size;
  *start_row = (size / numSections) * i;
  *end_row = (size / numSections) * (i + 1);
}

// helper function to initialize a size x size 2D array of int
int ** getIntMatrix(int size) {
  int ** counts = malloc(sizeof(int *) * size);
  if (counts == NULL) {
    printf("Error: failed malloc.  Exiting...\n");
    exit(1);
  }
  for (int i = 0; i < size; i++) {
    counts[i] = malloc(sizeof(int) * size);
    if (counts[i] == NULL) {
      printf("Error: failed malloc.  Exiting...\n");
      exit(1);
    }
  }
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      counts[i][j] = 0;
    }
  }
  return counts;
}

struct thread_data {
  int id;
  int size;
  int start_row;
  int end_row;
  int start_col;
  int end_col;
  float xmin;
  float xmax;
  float ymin;
  float ymax;
  int maxIterations_red;
  int maxIterations_green;
  int maxIterations_blue;
  struct ppm_pixel ** pixels;
  int ** in_set_red;
  int ** in_set_green;
  int ** in_set_blue;
  int ** counts_red;
  int ** counts_green;
  int ** counts_blue;
  int * max_count_red;  // largest count at a coordinate from counts
  int * max_count_green;
  int * max_count_blue;
};

void * thread_function(void * args) {
  struct thread_data * data = (struct thread_data *) args;
  printf("Thread %d) sub-image block: cols (%d, %d) to rows (%d, %d)\n",
      data->id, data->start_col, data->end_col, data->start_row, data->end_row);

  // compute membership for each color channel separately
  computeMembership(data->size, data->start_row, data->end_row, data->start_col,
      data->end_col, data->xmin, data->xmax, data->ymin, data->ymax,
      data->maxIterations_red, data->in_set_red);
  computeMembership(data->size, data->start_row, data->end_row, data->start_col,
      data->end_col, data->xmin, data->xmax, data->ymin, data->ymax,
      data->maxIterations_green, data->in_set_green);
  computeMembership(data->size, data->start_row, data->end_row, data->start_col,
      data->end_col, data->xmin, data->xmax, data->ymin, data->ymax,
      data->maxIterations_blue, data->in_set_blue);

  // compute counts for each color channel separately
  computeCounts(data->size, data->start_row, data->end_row, data->start_col,
      data->end_col, data->xmin, data->xmax, data->ymin, data->ymax,
      data->max_count_red, data->in_set_red, data->counts_red);
  computeCounts(data->size, data->start_row, data->end_row, data->start_col,
      data->end_col, data->xmin, data->xmax, data->ymin, data->ymax,
      data->max_count_green, data->in_set_green, data->counts_green);
  computeCounts(data->size, data->start_row, data->end_row, data->start_col,
      data->end_col, data->xmin, data->xmax, data->ymin, data->ymax,
      data->max_count_blue, data->in_set_blue, data->counts_blue);

  // wait for all threads before computing colors
  pthread_barrier_wait(&barrier);

  // compute colors for each color channel separately
  computeColors(data->start_row, data->end_row, data->start_col, data->end_col,
      data->max_count_red, data->counts_red, data->pixels, 1);
  computeColors(data->start_row, data->end_row, data->start_col, data->end_col,
      data->max_count_green, data->counts_green, data->pixels, 2);
  computeColors(data->start_row, data->end_row, data->start_col, data->end_col,
      data->max_count_blue, data->counts_blue, data->pixels, 3);

  printf("Thread %d) finished\n", data->id);
  return (void *) NULL;
}

//===========================================================//
//===========================================================//
//===========================================================//

int main(int argc, char* argv[]) {
  int size = 480;
  float xmin = -2.0;
  float xmax = 0.47;
  float ymin = -1.12;
  float ymax = 1.12;
  // each color channel has a different max iterations
  int maxIterations_red = 5000;
  int maxIterations_green = 500;
  int maxIterations_blue = 50;
  int numProcesses = 4;
  int start_col, end_col, start_row, end_row;
  struct ppm_pixel ** pixels = NULL;
  // whether a point is in the mandelbrot set for a color's maxIterations
  int ** in_set_red = NULL;
  int ** in_set_green = NULL;
  int ** in_set_blue = NULL;
  // each color channel keeps its own counts
  int ** counts_red = NULL;
  int ** counts_green = NULL;
  int ** counts_blue = NULL;
  double timer;
  struct timeval tstart, tend;
  char new_file[100];
  pthread_t * threads;
  struct thread_data * data;
  int ret1, ret2;  // for error checking
  // each color channel has its own max count
  int max_count_red = 0;
  int max_count_green = 0;
  int max_count_blue = 0;

  int opt;
  while ((opt = getopt(argc, argv, ":s:l:r:t:b:p:R:G:B:")) != -1) {
    switch (opt) {
      case 's': size = atoi(optarg); break;
      case 'l': xmin = atof(optarg); break;
      case 'r': xmax = atof(optarg); break;
      case 't': ymax = atof(optarg); break;
      case 'b': ymin = atof(optarg); break;
      case 'p': numProcesses = atoi(optarg); break;
      case 'R': maxIterations_red = atoi(optarg); break;
      case 'G': maxIterations_green = atoi(optarg); break;
      case 'B': maxIterations_blue = atoi(optarg); break;
      case '?': printf("usage: %s -s <size> -l <xmin> -r <xmax> "
        "-b <ymin> -t <ymax> -p <numProcesses> -R <maxIterationsRed> "
        "-G <maxIterationsGreen> -B <maxIterationsBlue>\n", argv[0]); break;
    }
  }
  printf("Generating color buddhabrot with size %dx%d\n", size, size);
  printf("  Num processes = %d\n", numProcesses);
  printf("  X range = [%.4f,%.4f]\n", xmin, xmax);
  printf("  Y range = [%.4f,%.4f]\n", ymin, ymax);

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

  // allocate memory to track set membership
  in_set_red = getIntMatrix(size);
  in_set_green = getIntMatrix(size);
  in_set_blue = getIntMatrix(size);

  // allocate memory for counts
  counts_red = getIntMatrix(size);
  counts_green = getIntMatrix(size);
  counts_blue = getIntMatrix(size);

  // initialize barrier and mutex
  ret1 = pthread_barrier_init(&barrier, NULL, numProcesses);
  if (ret1) {
    printf("ERROR: pthread_barrier_init failed\n");
    exit(0);
  }
  ret2 = pthread_mutex_init(&mutex, NULL);
  if (ret2) {
    printf("ERROR: pthread_mutex_init failed\n");
    exit(0);
  }

  gettimeofday(&tstart, NULL);

  // compute image
  for (int i = 0; i < numProcesses; i++) {
    getCoordinates(i, numProcesses, size, &start_col, &end_col,
        &start_row, &end_row);
    data[i].id = i;
    data[i].size = size;
    data[i].start_row = start_row;
    data[i].end_row = end_row;
    data[i].start_col = start_col;
    data[i].end_col = end_col;
    data[i].xmin = xmin;
    data[i].xmax = xmax;
    data[i].ymin = ymin;
    data[i].ymax = ymax;
    data[i].maxIterations_red = maxIterations_red;
    data[i].maxIterations_green = maxIterations_green;
    data[i].maxIterations_blue = maxIterations_blue;
    data[i].pixels = pixels;
    data[i].in_set_red = in_set_red;
    data[i].in_set_green = in_set_green;
    data[i].in_set_blue = in_set_blue;
    data[i].counts_red = counts_red;
    data[i].counts_green = counts_green;
    data[i].counts_blue = counts_blue;
    data[i].max_count_red = &max_count_red;
    data[i].max_count_green = &max_count_green;
    data[i].max_count_blue = &max_count_blue;
    // create threads
    pthread_create(&threads[i], NULL, thread_function, (void *) &data[i]);
  }

  // join threads
  for (int i = 0; i < numProcesses; i++) {
    pthread_join(threads[i], NULL);
  }

  gettimeofday(&tend, NULL);
  timer = tend.tv_sec - tstart.tv_sec + (tend.tv_usec - tstart.tv_usec)/1.e6;
  printf("Computed color buddhabrot set (%dx%d) in %g seconds\n", size, size, timer);

  // write to file
  new_file[0] = '\0';
  sprintf(new_file, "color_buddhabrot-%d-%d-%d-%d-%lu.ppm", size, maxIterations_red,
      maxIterations_green, maxIterations_blue, time(0));
  printf("Writing file %s\n", new_file);
  write_ppm(new_file, pixels, size, size);

  // free allocated array memory
  for (int i = 0; i < size; i++) {
    free(pixels[i]);
    pixels[i] = NULL;
    free(in_set_red[i]);
    in_set_red[i] = NULL;
    free(in_set_green[i]);
    in_set_green[i] = NULL;
    free(in_set_blue[i]);
    in_set_blue[i] = NULL;
    free(counts_red[i]);
    counts_red[i] = NULL;
    free(counts_green[i]);
    counts_green[i] = NULL;
    free(counts_blue[i]);
    counts_blue[i] = NULL;
  }
  free(pixels);
  pixels = NULL;
  free(in_set_red);
  in_set_red = NULL;
  free(in_set_green);
  in_set_green = NULL;
  free(in_set_blue);
  in_set_blue = NULL;
  free(counts_red);
  counts_red = NULL;
  free(counts_green);
  counts_green = NULL;
  free(counts_blue);
  counts_blue = NULL;
  free(threads);
  threads = NULL;
  free(data);
  data = NULL;

  // destroy barrier and mutex
  pthread_barrier_destroy(&barrier);
  pthread_mutex_destroy(&mutex);

  return 0;
}
