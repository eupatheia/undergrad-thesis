#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <assert.h>
#include <time.h>
#include <pthread.h>
#include <sys/time.h>
#include "read_ppm.h"
#include "mandelbrot_functions.h"

void getCoordinates(int i, int size, int * start_col, int * end_col,
    int * start_row, int * end_row);

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
  struct ppm_pixel * palette;
};

// function run by each thread to compute mandelbrot
void * thread_function(void * args) {
  struct thread_data * data = (struct thread_data *) args;
  printf("Thread %d) sub-image block: cols (%d, %d) to rows (%d, %d)\n",
      data->id, data->start_col, data->end_col, data->start_row, data->end_row);
  computeMandelbrot(data->size, data->start_row, data->end_row, data->start_col,
      data->end_col, data->xmin, data->xmax, data->ymin, data->ymax,
      data->maxIterations, data->pixels, data->palette);
  printf("Thread %d) finished\n", data->id);
  return (void *) NULL;
}

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
  struct ppm_pixel * palette = NULL;
  double timer;
  struct timeval tstart, tend;
  char new_file[100];
  pthread_t threads[4];
  struct thread_data data[4];

  // process command line arguments
  int opt;
  while ((opt = getopt(argc, argv, ":s:l:r:t:b:i:")) != -1) {
    switch (opt) {
      case 's': size = atoi(optarg); break;
      case 'l': xmin = atof(optarg); break;
      case 'r': xmax = atof(optarg); break;
      case 't': ymax = atof(optarg); break;
      case 'b': ymin = atof(optarg); break;
      case 'i': maxIterations = atoi(optarg); break;
      case '?': printf("usage: %s -s <size> -l <xmin> -r <xmax> "
        "-b <ymin> -t <ymax> -i <maxIterations>\n", argv[0]); break;
    }
  }
  printf("Generating mandelbrot with size %dx%d\n", size, size);
  printf("  Num processes = %d\n", numProcesses);
  printf("  X range = [%.4f,%.4f]\n", xmin, xmax);
  printf("  Y range = [%.4f,%.4f]\n", ymin, ymax);

  // attempt to allocate a 2D array of pixels
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

  // generate palette
  palette = malloc(sizeof(struct ppm_pixel) * maxIterations);
  if (palette == NULL) {
    printf("Error: failed malloc.  Exiting...\n");
    exit(1);
  }
  generatePalette(palette, maxIterations);

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
    data[i].palette = palette;
    // create threads
    pthread_create(&threads[i], NULL, thread_function, (void *) &data[i]);
  }

  // join threads
  for (int i = 0; i < numProcesses; i++) {
    pthread_join(threads[i], NULL);
  }

  gettimeofday(&tend, NULL);
  timer = tend.tv_sec - tstart.tv_sec + (tend.tv_usec - tstart.tv_usec)/1.e6;
  printf("Computed mandelbrot set (%dx%d) in %g seconds\n", size, size, timer);

  // write to file
  new_file[0] = '\0';
  sprintf(new_file, "mandelbrot_S%d_I%d_L%.3f_R%.3f_B%.3f_T%.3f_%lu.ppm",
      size, maxIterations, xmin, xmax, ymin, ymax, time(0));
  printf("Writing file %s\n", new_file);
  write_ppm(new_file, pixels, size, size);

  // free allocated array memory
  for (int i = 0; i < size; i++) {
    free(pixels[i]);
    pixels[i] = NULL;
  }
  free(pixels);
  pixels = NULL;
  free(palette);
  palette = NULL;

  return 0;
}

// assign quadrant coordinates
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
