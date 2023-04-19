#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <assert.h>
#include <time.h>
#include <pthread.h>
#include <sys/time.h>
#include <math.h>
#include "read_ppm.h"

// Step 1: Determine mandelbrot set membership
void computeMembership(int size, int start_row, int end_row, int start_col,
    int end_col, float xmin, float xmax, float ymin, float ymax, float cx,
    float cy, int maxIterations, struct ppm_pixel ** pixels) {

  float xfrac, yfrac, x, y, xtmp;
  int iter;

  for (int i = start_row; i < end_row; i++) {
    for (int j = start_col; j < end_col; j++) {
      // translate coordinates to center mandelbrot object
      xfrac = (float) j / size;
      yfrac = (float) i / size;
      // this coordinate on the "graph" will be the seed
      x = xmin + xfrac * (xmax - xmin);
      y = ymin + yfrac * (ymax - ymin);

      iter = 0;
      // square both to avoid using sqrt
      float cutoff = fmax(cx*cx + cy*cy, 2*2);
      while (iter < maxIterations && x*x + y*y < cutoff) {
        xtmp = x*x - y*y + cx;
        y = 2*x*y + cy;
        x = xtmp;
        iter++;
      }

      if (iter < maxIterations) {
        // escaped
        int section = maxIterations / 5;
        if (iter < section) {
          pixels[i][j].red = 255;
          pixels[i][j].green = round(((float) (iter % section) / section) * 255);
          pixels[i][j].blue = 0;
        } else if (iter < section * 2) {
          pixels[i][j].red = 255 - round(((float) (iter % section) / section) * 255);
          pixels[i][j].green = 255;
          pixels[i][j].blue = 0;
        } else if (iter < section * 3) {
          pixels[i][j].red = 0;
          pixels[i][j].green = 255;
          pixels[i][j].blue = round(((float) (iter % section) / section) * 255);
        } else if (iter < section * 4) {
          pixels[i][j].red = 0;
          pixels[i][j].green = 255 - round(((float) (iter % section) / section) * 255);
          pixels[i][j].blue = 255;
        } else {
          pixels[i][j].red = round(((float) (iter % section) / section) * 255);
          pixels[i][j].green = 0;
          pixels[i][j].blue = 255;
        }
      } else {
        // did not escape, color black
        pixels[i][j].red = 0;
        pixels[i][j].green = 0;
        pixels[i][j].blue = 0;
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
  int size;
  int start_row;
  int end_row;
  int start_col;
  int end_col;
  float xmin;
  float xmax;
  float ymin;
  float ymax;
  float cx;
  float cy;
  struct ppm_pixel ** pixels;
  int maxIterations;
};

void * thread_function(void * args) {
  struct thread_data * data = (struct thread_data *) args;
  printf("Thread %d) sub-image block: cols (%d, %d) to rows (%d, %d)\n",
      data->id, data->start_col, data->end_col, data->start_row, data->end_row);

  computeMembership(data->size, data->start_row, data->end_row, data->start_col,
      data->end_col, data->xmin, data->xmax, data->ymin, data->ymax, data->cx,
      data->cy, data->maxIterations, data->pixels);

  printf("Thread %d) finished\n", data->id);
  return (void *) NULL;
}

//===========================================================//
//===========================================================//
//===========================================================//

int main(int argc, char* argv[]) {
  int size = 480;
  float xmin = -2.0;
  float xmax = 2.0;
  float ymin = -2.0;
  float ymax = 2.0;
  float cx = -1.0;  // real component of parameter c
  float cy = 0.0;   // imaginary component of parameter c
  int numProcesses = 4;
  int maxIterations = 60;
  int start_col, end_col, start_row, end_row;
  struct ppm_pixel ** pixels = NULL;
  double timer;
  struct timeval tstart, tend;
  char new_file[100];
  pthread_t * threads;
  struct thread_data * data;

  int opt;
  while ((opt = getopt(argc, argv, ":s:l:r:t:b:p:x:y:i:")) != -1) {
    switch (opt) {
      case 's': size = atoi(optarg); break;
      case 'l': xmin = atof(optarg); break;
      case 'r': xmax = atof(optarg); break;
      // switch t and b, else -y axis points up
      case 't': ymin = atof(optarg); break;
      case 'b': ymax = atof(optarg); break;
      case 'p': numProcesses = atoi(optarg); break;
      case 'x': cx = atof(optarg); break;
      case 'y': cy = atof(optarg); break;
      case 'i': maxIterations = atoi(optarg); break;
      case '?': printf("usage: %s -s <size> -l <xmin> -r <xmax> "
        "-b <ymin> -t <ymax> -p <numProcesses> -x <c-real-comp> "
        "-y <c-imaginary-comp> -i <maxIterations>\n", argv[0]); break;
    }
  }
  printf("Generating quadratic Julia set with size %dx%d\n", size, size);
  printf("  maxIterations = %d\n", maxIterations);
  printf("  c = %.3f + %.3fi\n", cx, cy);
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
    data[i].cx = cx;
    data[i].cy = cy;
    data[i].pixels = pixels;
    data[i].maxIterations = maxIterations;
    // create threads
    pthread_create(&threads[i], NULL, thread_function, (void *) &data[i]);
  }

  // join threads
  for (int i = 0; i < numProcesses; i++) {
    pthread_join(threads[i], NULL);
  }

  gettimeofday(&tend, NULL);
  timer = tend.tv_sec - tstart.tv_sec + (tend.tv_usec - tstart.tv_usec)/1.e6;
  printf("Computed Julia set (%dx%d) in %g seconds\n", size, size, timer);

  // write to file
  new_file[0] = '\0';
  sprintf(new_file, "julia_S%d_CX%.3f_CY%.3f_I%d_L%.3f_R%.3f_B%.3f_T%.3f_%lu.ppm",
      size, cx, cy, maxIterations, xmin, xmax, ymin, ymax, time(0));
  printf("Writing file %s\n", new_file);
  write_ppm(new_file, pixels, size, size);

  // free allocated array memory
  for (int i = 0; i < size; i++) {
    free(pixels[i]);
    pixels[i] = NULL;
  }
  free(pixels);
  pixels = NULL;
  free(threads);
  threads = NULL;
  free(data);
  data = NULL;

  return 0;
}
