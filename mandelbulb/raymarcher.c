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
float stepSize = 0.1;
float softness = 8;  // softness factor for shadows

struct lighting {
  struct vec ambient;
  struct vec diffuse;
  struct vec specular;
};

// calculate SDF (signed distance function) from a position pos
// to a sphere defined at center center and radius radius 
float sdSphere(struct vec pos, struct vec center, float radius) {
  return distance(pos, center) - radius;
}

// calculate SDF of a position pos from a plane at height h with normal n
float sdPlane(struct vec pos, struct vec norm, float height) {
  norm = normalize(norm);
  return dot(pos, norm) + height;
}

// calculate min SDF to any object in scene
float sdScene(struct vec pos) {
  float sdf;
  struct vec c1 = {6, 0, 0, 1.0};
  float r1 = 3.0;
  struct vec c2 = {0, 3, 0, 1.0};
  float r2 = 2.0;
  sdf = fmin(sdSphere(pos, c1, r1), sdSphere(pos, c2, r2));
  struct vec planeNorm = {0, 1, 0, 0};
  sdf = fmin(sdf, sdPlane(pos, planeNorm, 10));
  return sdf;
}

// calculate the normal at some position pos on a surface
// by calculating the gradient in every direction (how does a small change
// by epsilon along one axis change the SDF?)
struct vec calcNormal(struct vec pos, float epsilon) {
  struct vec x1 = {pos.x - epsilon, pos.y, pos.z, 1.0};
  struct vec x2 = {pos.x + epsilon, pos.y, pos.z, 1.0};
  struct vec y1 = {pos.x, pos.y - epsilon, pos.z, 1.0};
  struct vec y2 = {pos.x, pos.y + epsilon, pos.z, 1.0};
  struct vec z1 = {pos.x, pos.y, pos.z - epsilon, 1.0};
  struct vec z2 = {pos.x, pos.y, pos.z + epsilon, 1.0};
  float deltaX = sdScene(x2) - sdScene(x1);
  float deltaY = sdScene(y2) - sdScene(y1);
  float deltaZ = sdScene(z2) - sdScene(z1);
  struct vec norm = {deltaX, deltaY, deltaZ, 0.0};
  return normalize(norm);  // return unit vector
}

// calculate phong shading colors
struct lighting phongShading(struct vec pos, struct vec norm, struct vec La,
    struct vec Ld, struct vec Ls, struct vec Ka, struct vec Kd,
    struct vec Ks, float shininess) {
  // normalized direction to light source
  struct vec s = normalize(vSub(lightPos, pos));
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
  struct lighting shades = {ambient, diffuse, specular};
  return shades;
}

// BUGGY(?) distance-aided raymarch:
// given a normalized ray with tail at pos,
// performs raymarch and returns 1 if hit, else 0,
// and returns position of surface in hitPos if hit, and
//             shadow brightness in [0, 1] (0 = dark, 1 bright);
// see https://iquilezles.org/articles/rmshadows/ for soft shadows
int DAraymarch(struct vec pos, struct vec ray, struct vec * hitPos,
    float * shadowColor) {
  *shadowColor = 1;  // start with full brightness
  float totalDistance = 0.0f;
  while (totalDistance < maxMarchDistance) {
    struct vec currPos = vAdd(pos, scale(ray, totalDistance));
    float minSDF = fabsf(sdScene(currPos));
    if (minSDF < hitRange) {
      // close enough, call a hit
      hitPos->x = currPos.x;
      hitPos->y = currPos.y;
      hitPos->z = currPos.z;
      hitPos->a = currPos.a;
      *shadowColor = 0.1;
      return 1;
    }
    // clamp between 0.1 and 1
    *shadowColor = fmax(fmin(*shadowColor, softness * minSDF / totalDistance), 0.1);
    // printf("%.3f = %.3f\n", minSDF, *shadowColor);
    totalDistance += minSDF;
  }
  // never hit any surface in scene
  return 0;
}

// REGULAR equal-step raymarch:
// given a normalized ray with tail at pos,
// performs raymarch and returns 1 if hit, else 0,
// and returns position of surface in hitPos if hit
int STEPraymarch(struct vec pos, struct vec ray, struct vec * hitPos,
    float * shadowColor) {
  *shadowColor = 1;  // start with full brightness
  for (float i = 0; i < maxMarchDistance; i += stepSize) {
    struct vec currPos = vAdd(pos, scale(ray, i));
    float minSDF = fabsf(sdScene(currPos));
    if (minSDF <= (stepSize / 2.0)) {
      // close enough, call a hit
      hitPos->x = currPos.x;
      hitPos->y = currPos.y;
      hitPos->z = currPos.z;
      hitPos->a = currPos.a;
      *shadowColor = 0;
      return 1;
    }
    *shadowColor = fmin(*shadowColor, softness * minSDF / i);
  }
  // never hit any surface in scene
  return 0;
}

// computes the color for a specific position and ray direction
struct ppm_pixel computeColor(struct vec pos, struct vec ray) {
  float shadowColor;
  struct ppm_pixel color;
  struct vec hitPos = {0, 0, 0, 0};
  // first compute primary ray
  int hit = DAraymarch(pos, ray, &hitPos, &shadowColor);
  if (hit == 1) {
    // calculate phong shading color first
    struct vec norm = calcNormal(hitPos, 0.0006);
    struct vec lightColor = {1.0, 1.0, 1.0, 0.0};
    // struct vec Ka = {0.24725f, 0.2245f, 0.0645f, 0.0f};
    // struct vec Kd = {0.34615f, 0.3143f, 0.0903f, 0.0f};
    // struct vec Ks = {0.797357f, 0.72399f, 0.20801f, 0.0f};
    struct vec Ka = {0.1, 0.1, 0.1, 0};
    struct vec Kd = {0.5, 0.5, 0.5, 0};
    struct vec Ks = {0.5, 0.5, 0.5, 0};
    struct lighting shades = phongShading(hitPos, norm, lightColor, lightColor,
        lightColor, Ka, Kd, Ks, 20.0f);
    // then calculate shadow ray
    // normalized direction to light source
    struct vec lightDir = normalize(vSub(lightPos, hitPos));
    // printf("%.3f %.3f %.3f\n", lightDir.x, lightDir.y, lightDir.z);
    struct vec shadowPos = {0, 0, 0, 0};
    // start slightly away from surface to avoid hitting same point again
    struct vec hitPosPlus = vAdd(hitPos, scale(lightDir, hitRange * 10));

    int hasShadow = DAraymarch(hitPosPlus, lightDir, &shadowPos, &shadowColor);
    struct vec finalColor;
    if (dot(norm, lightDir) <= 0) {
      // surface not facing light, add all components
      finalColor = vAdd(vAdd(shades.ambient, shades.diffuse),
          shades.specular);
    } else {
      // darken by shadow factor
      finalColor = scale(vAdd(vAdd(shades.ambient, shades.diffuse),
          shades.specular), shadowColor);
    }
    color.red = fmin(finalColor.x * 255, 255);
    color.green = fmin(finalColor.y * 255, 255);
    color.blue = fmin(finalColor.z * 255, 255);
    return color;
  } else {
    // primary ray never hit surface, color black
    color.red = 0;
    color.green = 0;
    color.blue = 0;
    return color;
  }
}

// render pixels from start rows/cols to end rows/cols
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
      pixels[i][j] = computeColor(camPos, ray);
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
  sprintf(new_file, "raymarcher_S%d_%lu.ppm", size, time(0));
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