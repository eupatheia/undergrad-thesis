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
    int * max_count, int ** counts, struct ppm_pixel ** pixels) {

  float gamma = 0.681;
  float factor = 1.0 / gamma;

  for (int i = start_row; i < end_row; i++) {
    for (int j = start_col; j < end_col; j++) {
      float value = 0;
      if (counts[i][j] > 0) {
	value = log(counts[i][j]) / log(*max_count);
        value = pow(value, factor);
      }
      pixels[i][j].red = value * 255;
      pixels[i][j].green = value * 255;
      pixels[i][j].blue = value * 255;
    }
  }
}

// helper function to assign quadrant coordinates
void getCoordinates(int i, int size, int * start_col, int * end_col,
    int * start_row, int * end_row) {
  if (i == 0) {
    // top left
    *start_col = 0;
    *end_col = size / 2;
    *start_row = 0;
    *end_row = size / 2;
  } else if (i == 1) {
    // top right
    *start_col = size / 2;
    *end_col = size;
    *start_row = 0;
    *end_row = size / 2;
  } else if (i == 2) {
    // bottom left
    *start_col = 0;
    *end_col = size / 2;
    *start_row = size / 2;
    *end_row = size;
  } else {
    // bottom right
    *start_col = size / 2;
    *end_col = size;
    *start_row = size / 2;
    *end_row = size;
  }
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
  int maxIterations;
  struct ppm_pixel ** pixels;
  int ** in_set;
  int ** counts;
  int * max_count;  // largest count at a coordinate from counts
};

void * thread_function(void * args) {
  struct thread_data * data = (struct thread_data *) args;
  printf("Thread %d) sub-image block: cols (%d, %d) to rows (%d, %d)\n",
      data->id, data->start_col, data->end_col, data->start_row, data->end_row);

  computeMembership(data->size, data->start_row, data->end_row, data->start_col,
      data->end_col, data->xmin, data->xmax, data->ymin, data->ymax,
      data->maxIterations, data->in_set);

  computeCounts(data->size, data->start_row, data->end_row, data->start_col,
      data->end_col, data->xmin, data->xmax, data->ymin, data->ymax,
      data->max_count, data->in_set, data->counts);

  // wait for all threads before computing colors
  pthread_barrier_wait(&barrier);

  computeColors(data->start_row, data->end_row, data->start_col, data->end_col,
      data->max_count, data->counts, data->pixels);

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
  int maxIterations = 1000;
  int numProcesses = 4;
  int start_col, end_col, start_row, end_row;
  struct ppm_pixel ** pixels = NULL;
  int ** in_set = NULL;  // whether a point is in the mandelbrot set
  int ** counts = NULL;
  double timer;
  struct timeval tstart, tend;
  char new_file[100];
  pthread_t threads[4];
  struct thread_data data[4];
  int ret1, ret2;  // for error checking
  int max_count = 0;

  int opt;
  while ((opt = getopt(argc, argv, ":s:l:r:t:b:p:")) != -1) {
    switch (opt) {
      case 's': size = atoi(optarg); break;
      case 'l': xmin = atof(optarg); break;
      case 'r': xmax = atof(optarg); break;
      case 't': ymax = atof(optarg); break;
      case 'b': ymin = atof(optarg); break;
      case '?': printf("usage: %s -s <size> -l <xmin> -r <xmax> "
        "-b <ymin> -t <ymax> -p <numProcesses>\n", argv[0]); break;
    }
  }
  printf("Generating buddhabrot with size %dx%d\n", size, size);
  printf("  Num processes = %d\n", numProcesses);
  printf("  X range = [%.4f,%.4f]\n", xmin, xmax);
  printf("  Y range = [%.4f,%.4f]\n", ymin, ymax);

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
  in_set = malloc(sizeof(int *) * size);
  if (in_set == NULL) {
    printf("Error: failed malloc.  Exiting...\n");
    exit(1);
  }
  for (int i = 0; i < size; i++) {
    in_set[i] = malloc(sizeof(int) * size);
    if (in_set[i] == NULL) {
      printf("Error: failed malloc.  Exiting...\n");
      exit(1);
    }
  }

  // allocate memory for counts
  counts = malloc(sizeof(int *) * size);
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

  // initialize barrier and mutex
  ret1 = pthread_barrier_init(&barrier, NULL, 4);
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
    getCoordinates(i, size, &start_col, &end_col, &start_row, &end_row);
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
    data[i].maxIterations = maxIterations;
    data[i].pixels = pixels;
    data[i].in_set = in_set;
    data[i].counts = counts;
    data[i].max_count = &max_count;
    // create threads
    pthread_create(&threads[i], NULL, thread_function, (void *) &data[i]);
  }

  // join threads
  for (int i = 0; i < numProcesses; i++) {
    pthread_join(threads[i], NULL);
  }

  gettimeofday(&tend, NULL);
  timer = tend.tv_sec - tstart.tv_sec + (tend.tv_usec - tstart.tv_usec)/1.e6;
  printf("Computed buddhabrot set (%dx%d) in %g seconds\n", size, size, timer);

  // write to file
  new_file[0] = '\0';
  sprintf(new_file, "buddhabrot-%d-%lu.ppm", size, time(0));
  printf("Writing file %s\n", new_file);
  write_ppm(new_file, pixels, size, size);

  // free allocated array memory
  for (int i = 0; i < size; i++) {
    free(pixels[i]);
    pixels[i] = NULL;
    free(in_set[i]);
    in_set[i] = NULL;
    free(counts[i]);
    counts[i] = NULL;
  }
  free(pixels);
  pixels = NULL;
  free(in_set);
  in_set = NULL;
  free(counts);
  counts = NULL;

  // destroy barrier and mutex
  pthread_barrier_destroy(&barrier);
  pthread_mutex_destroy(&mutex);

  return 0;
}
