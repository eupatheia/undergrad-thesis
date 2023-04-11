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
float xmin = -1.5;
float xmax = 1.5;
float ymin = -1.5;
float ymax = 1.5;
struct vec camPos = {0, 0, -5, 1.0};
struct vec lightPos = {-10, 10, -3, 1.0};
float hitRange = 0.001f;
float maxMarchDistance = 500.0f;
float stepSize = 0.1;
int maxSteps = 128;
int maxIterations = 4;

struct lighting {
  struct vec ambient;
  struct vec diffuse;
  struct vec specular;
};

// calculate distance to the closest intersection of ray with tail at pos and
// the sphere centered at center with radius radius;
// returns (near intersection, far intersection) in distance along ray
// from https://iquilezles.org/articles/intersectors/
struct vec sphereIntersect(struct vec pos, struct vec ray, struct vec center,
    float radius) {
  struct vec intersection;
  struct vec posDir = vSub(pos, center);
  // if b < 0, sphere is behind the ray
  float b = dot(posDir, ray);
  float c = dot(posDir, posDir) - (radius * radius);
  float h = (b * b) - c;
  if (h < 0.0) {
    // no intersection
    intersection.x = -1;
    intersection.y = -1;
  } else {
    h = sqrt(h);
    intersection.x = -b - h;
    intersection.y = -b + h;
  }
  return intersection;
}

// compute of next iteration of f, given previous v,
// adapted from https://iquilezles.org/articles/mandelbulb/:
// finds the eighth power of v, doubles the angle and adds back to origPos
struct vec f(struct vec origPos, struct vec v) {
  // cartesian to spherical coordinates
  float r = length(v);
  float theta = acos(v.y / r);
  float phi = atan2(v.x, v.z);

  // radius to 8th power
  r = pow(r, 8.0);
  // multiply all angles by 8
  theta = theta * 8.0;
  phi = phi * 8.0;

  // convert back to cartesian coordinates
  struct vec w;
  // using origPos as offset c?
  w.x = origPos.x + r * sin(theta) * sin(phi);
  w.y = origPos.y + r * cos(theta);
  w.z = origPos.z + r * sin(theta) * cos(phi);

  return w;
}

// optimized computation of next iteration of f, given previous v,
// adapted from https://iquilezles.org/articles/mandelbulb/:
// finds the eighth power of v, doubles the angle, and adds c
struct vec optf(struct vec origPos, struct vec v) {
  // replace pow(v, 8) by individual multiplications
  float x = v.x;
  float x2 = x * x;
  float x4 = x2 * x2;
  float y = v.y;
  float y2 = y * y;
  float y4 = y2 * y2;
  float z = v.z;
  float z2 = z * z;
  float z4 = z2 * z2;

  // replace trig functions
  float k3 = x2 + z2;
  float k2 = 1.0 / sqrt(k3 * k3 * k3 * k3 * k3 * k3 * k3);
  float k1 = x4 + y4 + z4 - (6.0 * y2 * z2) - (6.0 * x2 * y2) + (2.0 * z2 * x2);
  float k4 = x2 - y2 + z2;

  struct vec w;
  w.x = origPos.x + 64.0 * x * y * z * (x2 - z2) * k4 *
      (x4 - 6.0 * x2 * z2 + z4) * k1 * k2;
  w.y = origPos.y + (-16.0 * y2 * k3 * k4 * k4) + (k1 * k1);
  w.z = origPos.z + -8.0 * y * k4 * ((x4 * x4) - (28.0 * x4 * x2 * z2) +
      (70.0 * x4 * z4) - (28.0 * x2 * z2 * z4) + (z4 * z4)) * k1 * k2;

  return w;
}

// calculate SDF (signed distance function) from a position pos
// to the mandelbulb object, returning distance
// adapted from https://iquilezles.org/articles/mandelbulb/:
float sdScene(struct vec pos) {
  struct vec w = pos;
  float dw = 1.0;  // gradient of potential
  float m = length(w);  // modulus of a point (i.e. length, |w|)
  int iter = 0;
  while (iter < maxIterations && m < 16) {
    // compute gradient dw_{n+1} = 8 * |w_n|^7 * dw_n
    dw = 8.0 * pow(m, 7) * dw + 1.0;
    // get next iteration
    w = f(pos, w);
    m = length(w);
    iter++;
  }
  // calc distance d = |w| * log|w| / |dw|;
  // note this is modulus, not abs val;
  // multiply 0.25 to offset errors in sdf
  // (see https://iquilezles.org/articles/distancefractals/)
  return (m * log(m) / dw) * 0.25;
}

// calculate the normal at some position pos on a surface
// by calculating the gradient in every direction (how does a change
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

// distance-aided raymarch:
// given a normalized ray with tail at pos,
// performs raymarch and returns 1 if hit, else 0,
// and if hit, returns position of surface in hitPos
int DAraymarch(struct vec pos, struct vec ray, struct vec * hitPos) {
  // compute distance to bounding sphere to only search w/in bounds
  struct vec center = {0, 0, 0, 0};
  struct vec sphereDistance = sphereIntersect(pos, ray, center, 1.25f);
  if (sphereDistance.y < 0) {
    // no far intersection, never hit bounding sphere, never hit mandelbulb
    return 0;
  }
  float totalDistance = sphereDistance.x;  // near intersection with sphere
  int steps = 0;
  // stop when too many steps or went past far intersection
  while (steps < maxSteps && totalDistance <= sphereDistance.y) {
    struct vec currPos = vAdd(pos, scale(ray, totalDistance));
    float minSDF = sdScene(currPos);
    if (minSDF < hitRange) {
      // close enough, call a hit
      hitPos->x = currPos.x;
      hitPos->y = currPos.y;
      hitPos->z = currPos.z;
      hitPos->a = currPos.a;
      return 1;
    }
    totalDistance += minSDF;
    steps++;
  }
  // never hit mandelbulb surface in scene
  return 0;
}

// REGULAR equal-step raymarch:
// given a normalized ray with tail at pos,
// performs raymarch and returns 1 if hit, else 0,
// and if hit, returns position of surface in hitPos
int STEPraymarch(struct vec pos, struct vec ray, struct vec * hitPos) {
  float d = 0;
  for (float i = 0; i < maxMarchDistance; i += stepSize) {
    struct vec currPos = vAdd(pos, scale(ray, i));
    float d = sdScene(currPos);
    d = fabsf(d);
    if (d <= (stepSize / 2.0)) {
      // close enough, call a hit
      hitPos->x = currPos.x;
      hitPos->y = currPos.y;
      hitPos->z = currPos.z;
      hitPos->a = currPos.a;
      return 1;
    }
  }
  // never hit any surface in scene
  return 0;
}

// computes the color for a specific position and ray direction
struct ppm_pixel computeColor(struct vec pos, struct vec ray) {
  struct ppm_pixel color;
  struct vec hitPos = {0, 0, 0, 0};
  // first compute primary ray intersection with surface and normal
  int hit = DAraymarch(pos, ray, &hitPos);
  struct vec norm = calcNormal(hitPos, 0.0005);
  if (hit == 1) {
    // calculate phong shading color first
    struct vec lightColor = {1.0, 1.0, 1.0, 0.0};
    struct vec Ka = {0.24725f, 0.2245f, 0.0645f, 0.0f};
    struct vec Kd = {0.34615f, 0.3143f, 0.0903f, 0.0f};
    struct vec Ks = {0.797357f, 0.72399f, 0.20801f, 0.0f};
    // struct vec Ka = {0, 0, 0.2f, 0};
    // struct vec Kd = {0, 0.3, 0.7, 0};
    // struct vec Ks = {1, 1, 1, 0};
    struct lighting shades = phongShading(hitPos, norm, lightColor, lightColor,
        lightColor, Ka, Kd, Ks, 20.0f);

    // then calculate shadow ray
    // normalized direction to light source
    struct vec lightDir = normalize(vSub(lightPos, hitPos));
    struct vec shadowPos = {0, 0, 0, 0};
    // start slightly away from surface to avoid hitting same point again
    struct vec hitPosPlus = vAdd(hitPos, scale(lightDir, hitRange * 2));

    int hasShadow = DAraymarch(hitPosPlus, lightDir, &shadowPos);
    struct vec finalColor;
    if (hasShadow == 0) {
      // not in shadow, add all components
      finalColor = vAdd(vAdd(shades.ambient, shades.diffuse),
          shades.specular);
    } else {
      // add just ambient and diffuse components and darken by shadow factor
      finalColor = vAdd(shades.ambient, shades.diffuse);
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
  while ((opt = getopt(argc, argv, ":s:p:l:r:b:t:")) != -1) {
    switch (opt) {
      case 's': size = atoi(optarg); break;
      case 'p': numProcesses = atoi(optarg); break;
      case 'l': xmin = atoi(optarg); break;
      case 'r': xmax = atoi(optarg); break;
      case 'b': ymin = atoi(optarg); break;
      case 't': ymax = atoi(optarg); break;
      case '?': printf("usage: %s -s <size> -p <numProcesses> "
          "-l <xmin> -r <xmax> -b <ymin> -t <ymax>\n", argv[0]);
          break;
    }
  }
  printf("Generating image with size %dx%d\n", size, size);
  printf("  Num processes = %d\n", numProcesses);
  printf("  x = [%.3f,%.3f]\n", xmin, xmax);
  printf("  y = [%.3f,%.3f]\n", ymin, ymax);

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
  sprintf(new_file, "mandelbulb_S%d_L%.3f_R%.3f_B%.3f_T%.3f_%lu.ppm", size,
      xmin, xmax, ymin, ymax, time(0));
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