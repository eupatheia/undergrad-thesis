#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <assert.h>
#include <time.h>
#include <sys/wait.h>
#include <sys/time.h>
#include <sys/shm.h>
#include <sys/ipc.h>
#include "read_ppm.h"
#include "mandelbrot_functions.h"

void getCoordinates(int i, int size, int * start_col, int * end_col,
    int * start_row, int * end_row);

int main(int argc, char* argv[]) {
  int size = 480;
  float xmin = -2.0;
  float xmax = 0.47;
  float ymin = -1.12;
  float ymax = 1.12;
  int maxIterations = 1000;
  int numProcesses = 4;
  int shmid1, shmid2;
  int * shmids = NULL;
  int start_col, end_col, start_row, end_row;
  struct ppm_pixel ** pixels = NULL;
  struct ppm_pixel * palette = NULL;
  double timer;
  struct timeval tstart, tend;
  char new_file[100];

  // process command line arguments
  int opt;
  while ((opt = getopt(argc, argv, ":s:l:r:t:b:p:")) != -1) {
    switch (opt) {
      case 's': size = atoi(optarg); break;
      case 'l': xmin = atof(optarg); break;
      case 'r': xmax = atof(optarg); break;
      case 't': ymax = atof(optarg); break;
      case 'b': ymin = atof(optarg); break;
      case 'p': numProcesses = atof(optarg); break;
      case '?': printf("usage: %s -s <size> -l <xmin> -r <xmax> "
        "-b <ymin> -t <ymax> -p <numProcesses>\n", argv[0]); break;
    }
  }
  printf("Generating mandelbrot with size %dx%d\n", size, size);
  printf("  Num processes = %d\n", numProcesses);
  printf("  X range = [%.4f,%.4f]\n", xmin, xmax);
  printf("  Y range = [%.4f,%.4f]\n", ymin, ymax);

  // attempt to allocate shared memory
  shmid1 = shmget(IPC_PRIVATE, sizeof(struct ppm_pixel *) * size, 0644 | IPC_CREAT);
  if (shmid1 == -1) {
    perror("Error: cannot initialize outer shared memory\n");
    exit(1);
  }
  // attempt to attach shared memory
  pixels = shmat(shmid1, NULL, 0);
  if (pixels == (void *) -1) {
    perror("Error: cannot access outer shared memory\n");
    exit(1);
  }
  // allocate and attach the inner arrays
  shmids = malloc(sizeof(int) * size);
  if (shmids == NULL) {
    printf("Error: failed malloc.  Exiting...\n");
    exit(1);
  }
  for (int i = 0; i < size; i++) {
    shmids[i] = shmget(IPC_PRIVATE, sizeof(struct ppm_pixel) * size, 0644 | IPC_CREAT);
    if (shmids[i] == -1) {
      perror("Error: cannot initialize inner shared memory\n");
      exit(1);
    }
    pixels[i] = shmat(shmids[i], NULL, 0);
    if (pixels[i] == (void *) -1) {
      perror("Error: cannot access inner shared memory\n");
      exit(1);
    }
  }

  // generate palette
  shmid2 = shmget(IPC_PRIVATE, sizeof(struct ppm_pixel) * maxIterations, 0644 | IPC_CREAT);
  if (shmid2 == -1) {
    perror("Error: cannot initialize palette shared memory\n");
    exit(1);
  }
  palette = shmat(shmid2, NULL, 0);
  if (palette == (void *) -1) {
    perror("Error: cannot access palette shared memory\n");
    exit(1);
  }
  generatePalette(palette, maxIterations);

  gettimeofday(&tstart, NULL);

  // compute image
  for (int i = 0; i < numProcesses; i++) {
    int pid = fork();
    if (pid == 0) {
      getCoordinates(i, size, &start_col, &end_col, &start_row, &end_row);
      printf("%d) Sub-image block: cols (%d, %d) to rows (%d, %d)\n", getpid(),
          start_col, end_col, start_row, end_row);
      for (int i = start_row; i < end_row; i++) {
        for (int j = start_col; j < end_col; j++) {
          computeMandelbrot(size, i, j, xmin, xmax, ymin, ymax, maxIterations,
              pixels, palette);
        }
      }
      free(shmids);
      shmids = NULL;
      exit(0);
    } else {
      printf("Launched child process: %d\n", pid);
    }
  }

  // wait for all child processes to complete
  for (int i = 0; i < numProcesses; i++) {
    int status;
    int pid = wait(&status);
    printf("Child process complete: %d\n", pid);
  }

  gettimeofday(&tend, NULL);
  timer = tend.tv_sec - tstart.tv_sec + (tend.tv_usec - tstart.tv_usec)/1.e6;
  printf("Computed mandelbrot set (%dx%d) in %g seconds\n", size, size, timer);

  // write to file
  new_file[0] = '\0';
  sprintf(new_file, "multi-mandelbrot-%d-%lu.ppm", size, time(0));
  printf("Writing file %s\n", new_file);
  write_ppm(new_file, pixels, size, size);

  // detach and delete all shared memory
  for (int i = 0; i < size; i++) {
    if (shmdt(pixels[i]) == -1) {
      perror("Error: cannot detatch from inner shared memory\n");
      exit(1);
    }
    if (shmctl(shmids[i], IPC_RMID, 0) == -1) {
      perror("Error: cannot remove inner shared memory\n");
      exit(1);
    }
  }
  if (shmdt(pixels) == -1 || shmdt(palette)) {
    perror("Error: cannot detatch from outer shared memory\n");
    exit(1);
  }
  if (shmctl(shmid1, IPC_RMID, 0) == -1 || shmctl(shmid2, IPC_RMID, 0) == -1) {
    perror("Error: cannot remove outer shared memory\n");
    exit(1);
  }

  free(shmids);
  shmids = NULL;

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
