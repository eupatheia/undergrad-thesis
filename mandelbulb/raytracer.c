#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <pthread.h>
#include <sys/time.h>
#include <math.h>
#include "read_write.h"
#include "vec.h"

int size = 300;
int numProcesses = 4;
float xmin = -10.0;
float xmax = 10.0;
float ymin = -10.0;
float ymax = 10.0;
struct vec camPos = {0, 0, -10, 1.0};
struct vec lightPos = {-10, 10, -5, 1.0};
float hitRange = 0.001f;
float maxMarchDistance = 256.0f;

// calculate SDF (signed distance function) from a position pos
// to a sphere defined at center center and radius radius 
float sdSphere(struct vec pos, struct vec center, float radius) {
  return distance(pos, center) - radius;
}

// calculate min SDF to any object in scene
float sdScene(struct vec pos) {
  // just one sphere in scene for now
  struct vec c1 = {6, 0, 0, 1.0};
  float r1 = 3.0;
  struct vec c2 = {0, 3, 0, 1.0};
  float r2 = 2.0;
  return fmin(sdSphere(pos, c1, r1), sdSphere(pos, c2, r2));
}

// calculate the normal at some position pos on a surface
// by calculating the gradient in every direction (how does a small change
// along one axis change the SDF?)
struct vec calcNormal(struct vec pos) {
  struct vec x1 = {pos.x - 0.001, pos.y, pos.z, 1.0};
  struct vec x2 = {pos.x + 0.001, pos.y, pos.z, 1.0};
  struct vec y1 = {pos.x, pos.y - 0.001, pos.z, 1.0};
  struct vec y2 = {pos.x, pos.y + 0.001, pos.z, 1.0};
  struct vec z1 = {pos.x, pos.y, pos.z - 0.001, 1.0};
  struct vec z2 = {pos.x, pos.y, pos.z + 0.001, 1.0};
  float deltaX = sdScene(x2) - sdScene(x1);
  float deltaY = sdScene(y2) - sdScene(y1);
  float deltaZ = sdScene(z2) - sdScene(z1);
  struct vec norm = {deltaX, deltaY, deltaZ, 0.0};
  return normalize(norm);  // return unit vector
}

struct ppm_pixel raymarch(struct vec pos, struct vec ray);

// calculate phong shading color
struct vec phongShading(struct vec pos, struct vec norm, struct vec La,
    struct vec Ld, struct vec Ls, struct vec Ka, struct vec Kd,
    struct vec Ks, float shininess) {
  // normalized direction to light source
  struct vec s = normalize(vSub(lightPos, pos));
  // struct vec posPlusE = {pos.x + 0.000001, pos.y + 0.000001,
  //     pos.z + 0.000001, 0};
  // struct ppm_pixel shadowColor = raymarch(posPlusE, s);
  // normalized direction to viewer
  struct vec v = {camPos.x - pos.x, camPos.y - pos.y, camPos.z - pos.z, 0.0f};
  v = normalize(v);
  float sDotN = dot(s, norm);
  struct vec twoSDotN = {2.0 * sDotN, 2.0 * sDotN, 2.0 * sDotN, 0.0f};
  // reflection of s about norm
  struct vec r = vSub(vMul(twoSDotN, norm), s);
  struct vec ambient = vMul(Ka, La);
  struct vec diffuse = {fmax(sDotN, 0.0) * Ld.x * Kd.x,
      fmax(sDotN, 0.0) * Ld.y * Kd.y, fmax(sDotN, 0.0) * Ld.z * Kd.z, 0.0f};
  float shine = pow(fmax(dot(v, r), 0.0), shininess);
  struct vec specular = {Ks.x * Ls.x * shine, Ks.y * Ls.y * shine,
      Ks.z * Ls.z * shine, 0.0f};
  return vAdd(vAdd(ambient, diffuse), specular);
}

// given a normalized ray with tail at pos,
// performs raymarch and returns correct color
struct ppm_pixel raymarch(struct vec pos, struct vec ray) {
  struct ppm_pixel color;
  float totalDistance = 0.0f;
  while (totalDistance >= 0 && totalDistance < maxMarchDistance) {
    struct vec distanceVec = {totalDistance, totalDistance, totalDistance, 0.0};
    struct vec currPos = vAdd(pos, vMul(distanceVec, ray));
    float minSDF = sdScene(currPos);
    if (minSDF < hitRange) {
      // close enough, color surface
      struct vec norm = calcNormal(currPos);
      struct vec lightColor = {1.0, 1.0, 1.0, 0.0};
      // struct vec Ka = {0.24725f, 0.2245f, 0.0645f, 0.0f};
      // struct vec Kd = {0.34615f, 0.3143f, 0.0903f, 0.0f};
      // struct vec Ks = {0.797357f, 0.72399f, 0.20801f, 0.0f};
      struct vec Ka = {0, 0, 0.2f, 0};
      struct vec Kd = {0, 0.3, 0.7, 0};
      struct vec Ks = {1, 1, 1, 0};
      struct vec c = phongShading(currPos, norm, lightColor, lightColor,
          lightColor, Ka, Kd, Ks, 20.0f);
      color.red = fmin(c.x * 255, 255);
      color.green = fmin(c.y * 255, 255);
      color.blue = fmin(c.z * 255, 255);
      return color;
    }
    totalDistance += minSDF;
  }
  // never hit any surface in scene, color black
  color.red = 0;
  color.green = 0;
  color.blue = 0;
  return color;
}

// color pixels from start rows/cols to end rows/cols
void render(int start_row, int end_row, int start_col, int end_col,
    struct ppm_pixel ** pixels) {
  for (int i = start_row; i < end_row; i++) {
    for (int j = start_col; j < end_col; j++) {
      // convert row-col position to position in [xmin, xmax] and [ymin, ymax]
      float u = ((((float) j + 0.5) / size) * (xmax - xmin)) + xmin;
      float v = ((((float) i + 0.5) / size) * (ymax - ymin)) + ymin;
      v = -v;  // to make y increase upwards (otherwise image is flipped)
      struct vec ray = {u, v, -camPos.z, 0.0};
      ray = normalize(ray);  // direction as a unit vector
      pixels[i][j] = raymarch(camPos, ray);
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
};

void * thread_function(void * args) {
  struct thread_data * data = (struct thread_data *) args;
  printf("Thread %d) sub-image block: cols (%d, %d) to rows (%d, %d)\n",
      data->id, data->start_col, data->end_col, data->start_row, data->end_row);

  render(data->start_row, data->end_row, data->start_col, data->end_col,
      data->pixels);

  printf("Thread %d) finished\n", data->id);
  return (void *) NULL;
}

//===========================================================//
//===========================================================//
//===========================================================//

int main(int argc, char* argv[]) {
  int start_col, end_col, start_row, end_row;
  struct ppm_pixel ** pixels = NULL;
  double timer;
  struct timeval tstart, tend;
  char new_file[100];
  pthread_t * threads;
  struct thread_data * data;

  int opt;
  while ((opt = getopt(argc, argv, ":s:p:")) != -1) {
    switch (opt) {
      case 's': size = atoi(optarg); break;
      case 'p': numProcesses = atoi(optarg); break;
      case '?': printf("usage: %s -s <size> -p <numProcesses>\n", argv[0]);
          break;
    }
  }
  printf("Generating image with size %dx%d\n", size, size);
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
    data[i].start_row = start_row;
    data[i].end_row = end_row;
    data[i].start_col = start_col;
    data[i].end_col = end_col;
    data[i].pixels = pixels;
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
  sprintf(new_file, "raytracer_S%d_%lu.ppm", size, time(0));
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