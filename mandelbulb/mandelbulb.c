#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <pthread.h>
#include <sys/time.h>
#include <math.h>
#include "read_write.h"
#include "vec.h"
#include "mat.h"
#include "shaders.h"
#include "ray_functions.h"

// computes the color for a specific position (eyePos) and ray direction,
// using a specific shader
struct ppm_pixel computeColor(int shader, struct vec pos, struct vec ray,
    struct vec lightPos, int maxSteps, float hitRange, int maxIterations,
    float softness, float yratio) {
  struct ppm_pixel color;
  float shadowColor;
  struct vec hitPos = {0, 0, 0, 0};
  // first compute primary ray intersection with surface and normal
  int hit = DAraymarch(pos, ray, &hitPos, &shadowColor, maxSteps, hitRange,
      maxIterations, softness);
  // calculate surface normal at hitPos
  struct vec norm = calcNormal(hitPos, 0.0005, maxIterations);

  if (hit == 1) {
    // on the surface, pick between shaders
    if (shader == 1) {
      return phongShadowShader(pos, hitPos, norm, lightPos, hitRange, maxSteps,
          maxIterations, softness);
    } else if (shader == 2) {
      return reflectShader(ray, norm, hitPos);
    } else if (shader == 3) {
      return refractShader(ray, norm, hitPos, 0.9);
    } else if (shader == 4) {
      return dielectricShader(ray, norm, hitPos, 0.9);
    } else if (shader == 5) {
      struct vec ratios = {1.5, 2, 2.4, 1};
      return chromaticDispersionShader(ray, norm, hitPos, ratios);
    } else if (shader == 6) {
      return volumetricShader(pos, ray, hitPos, norm, yratio, maxIterations,
          0.9, maxSteps, hitRange, softness);
    } else {
      return normalShader(yratio, norm);  // default normal coloring
    }
  } else {
    // primary ray never hit surface, color background
    return backgroundColor(yratio);
  }
}

// render pixels from start rows/cols to end rows/cols
void render(int start_row, int end_row, int start_col, int end_col,
    struct ppm_pixel ** pixels, int size, float xmin, float xmax, float ymin,
    float ymax, struct mat VInv, struct vec camPos, struct vec lightPos,
    float hitRange, int maxSteps, int maxIterations, int AAWidth,
    float softness, int shader) {
  for (int i = start_row; i < end_row; i++) {
    for (int j = start_col; j < end_col; j++) {
      // multiple samples inside a pixel for anti-aliasing
      float inc = 1.0 / (AAWidth + 1);  // increments between samples in pixel
      int reds = 0;
      int greens = 0;
      int blues = 0;
      for (int m = 1; m <= AAWidth; m++) {
        for (int n = 1; n <= AAWidth; n++) {
          // convert row-col to position in [xmin, xmax] and [ymin, ymax]
          float u = ((((float) j + (n * inc)) / size) * (xmax - xmin)) + xmin;
          // to make y increase upwards (otherwise image is flipped)
          float v = ((((float) (size - i - 1) + (m * inc)) / size) *
                               (ymax - ymin)) + ymin;
          // ray position in eye coords (i.e. WRT camera)
          struct vec ray = {u, v, -5, 0};
          ray = trans(VInv, ray);
          ray = normalize(ray);
          if (i == 150 && j == 150) {
            printMat(VInv);
            printVec(ray);
          }
          struct ppm_pixel color = computeColor(shader, camPos, ray, lightPos,
              maxSteps, hitRange, maxIterations, softness, (float) i / size);
          reds += color.red;
          greens += color.green;
          blues += color.blue;
        }
      }
      // average color of all subsamples
      pixels[i][j].red = (float) reds / (AAWidth * AAWidth);
      pixels[i][j].green = (float) greens / (AAWidth * AAWidth);
      pixels[i][j].blue = (float) blues / (AAWidth * AAWidth);
    }
  }
}

// helper function to divide plane into sections and assign coordinates,
// where each section is like a "stripe" across the plane,
// and i is a section number from 0 to numSections - 1;
// a section is defined from start (inclusive) to end (exclusive)
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
  int size;
  float xmin;
  float xmax;
  float ymin;
  float ymax;
  struct mat VInv;
  struct vec camPos;
  struct vec lightPos;
  float hitRange;
  float maxMarchDistance;
  int maxSteps;
  float stepSize;
  int maxIterations;
  int AAWidth;
  float softness;
  int shader;
};

void * thread_function(void * args) {
  struct thread_data * data = (struct thread_data *) args;
  printf("Thread %d) sub-image block: cols (%d, %d) to rows (%d, %d)\n",
      data->id, data->start_col, data->end_col, data->start_row, data->end_row);

  render(data->start_row, data->end_row, data->start_col, data->end_col,
      data->pixels, data->size, data->xmin, data->xmax, data->ymin,
      data->ymax, data->VInv, data->camPos, data->lightPos,
      data->hitRange, data->maxSteps, data->maxIterations, data->AAWidth,
      data->softness, data->shader);

  printf("Thread %d) finished\n", data->id);
  return (void *) NULL;
}

//===========================================================//
//===========================================================//
//===========================================================//

int main(int argc, char* argv[]) {
  int size = 300;
  int numProcesses = 4;
  float xmin = -1.5;
  float xmax = 1.5;
  float ymin = -1.5;
  float ymax = 1.5;
  struct vec camPos = {0, 0, 5, 1};
  struct vec lookPos = {0, 0, 0, 1};  // look at origin
  struct vec upDir = {0, 1, 0, 0};
  struct vec lightPos = {-10, 10, 5, 1};
  float hitRange = 0.001f;
  float maxMarchDistance = 500.0f;
  float stepSize = 0.1;
  int maxSteps = 128;
  int maxIterations = 4;
  struct mat VInv;  // inverse view matrix
  int start_col, end_col, start_row, end_row;
  struct ppm_pixel ** pixels = NULL;
  double timer;
  struct timeval tstart, tend;
  char new_file[100];
  pthread_t * threads;
  struct thread_data * data;
  int AAWidth = 1;
  float softness = 8;
  int shader = 0;

  int opt;
  while ((opt = getopt(argc, argv, ":s:p:a:c:")) != -1) {
    switch (opt) {
      case 's': size = atoi(optarg); break;
      case 'p': numProcesses = atoi(optarg); break;
      case 'a': AAWidth = atoi(optarg); break;
      case 'c': shader = atoi(optarg); break;
      case '?': printf("usage: %s -s <size> -p <num processes> "
          "-a <antialiasingFilterWidth> -c <colorSetting>\n"
          "COLOR/SHADER SETTINGS:\n"
          "  1 = phong\n"
          "  2 = total reflection\n"
          "  3 = total refraction\n"
          "  4 = dielectric\n"
          "  5 = chromatic dispersion\n"
          "  6 = volumetric\n"
          "  (other) = normal\n", argv[0]);
          exit(0);
    }
  }

  // int radius = 5;
  // int totalImages = 36;
  // int degChange = 360 / totalImages;
  // for (int i = 0; i < 360; i += degChange) {
  //   float theta = i * (M_PI / 180.0); // in radians
  //   float camZ = radius * cos(theta);
  //   float camX = radius * sin(theta);
  //   camPos.x = camX;
  //   camPos.z = camZ;

    printf("Generating image with size %dx%d\n", size, size);
    printf("  Num processes = %d\n", numProcesses);
    printf("  Anti-aliasing filter width = %d\n", AAWidth);
    printf("  Shader setting = %d\n", shader);
    printf("  camPos = ");
    printVec(camPos);

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

    // calculate inverse matrices
    VInv = getInvViewMat(camPos, lookPos, upDir);

    // compute image
    for (int i = 0; i < numProcesses; i++) {
      getCoordinates(i, numProcesses, size, &start_col, &end_col,
          &start_row, &end_row);
      data[i].id = i;
      data[i].start_row = start_row;
      data[i].end_row = end_row;
      data[i].start_col = start_col;
      data[i].end_col = end_col;
      data[i].pixels = pixels;
      data[i].size = size;
      data[i].xmin = xmin;
      data[i].xmax = xmax;
      data[i].ymin = ymin;
      data[i].ymax = ymax;
      data[i].VInv = VInv;
      data[i].camPos = camPos;
      data[i].lightPos = lightPos;
      data[i].hitRange = hitRange;
      data[i].maxMarchDistance = maxMarchDistance;
      data[i].maxSteps = maxSteps;
      data[i].stepSize = stepSize;
      data[i].maxIterations = maxIterations;
      data[i].AAWidth = AAWidth;
      data[i].softness = softness;
      data[i].shader = shader;
      // create threads
      pthread_create(&threads[i], NULL, thread_function, (void *) &data[i]);
    }

    // join threads
    for (int i = 0; i < numProcesses; i++) {
      pthread_join(threads[i], NULL);
    }

    gettimeofday(&tend, NULL);
    timer = tend.tv_sec - tstart.tv_sec + (tend.tv_usec - tstart.tv_usec)/1.e6;
    printf("Computed image (%dx%d) in %g seconds\n", size, size, timer);

    // write to file
    new_file[0] = '\0';
    sprintf(new_file, "mandelbulb_S%d_A%d_C%d_%lu.ppm", size, AAWidth, shader,
        time(0));
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
  // }

  return 0;
}