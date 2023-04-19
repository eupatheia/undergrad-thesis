#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "vec.h"
#include "mat.h"

void printMat(const struct mat m) {
  printf("  |%.3f\t%.3f\t%.3f\t%.3f\t|\n", m.a.x, m.b.x, m.c.x, m.d.x);
  printf("  |%.3f\t%.3f\t%.3f\t%.3f\t|\n", m.a.y, m.b.y, m.c.y, m.d.y);
  printf("  |%.3f\t%.3f\t%.3f\t%.3f\t|\n", m.a.z, m.b.z, m.c.z, m.d.z);
  printf("  |%.3f\t%.3f\t%.3f\t%.3f\t|\n\n", m.a.a, m.b.a, m.c.a, m.d.a);
}

// multiply two 4x4 matrices to get another 4x4 matrix
struct mat matMul(const struct mat m, const struct mat n) {
  struct mat p;
  struct vec r1 = {m.a.x, m.b.x, m.c.x, m.d.x};
  struct vec r2 = {m.a.y, m.b.y, m.c.y, m.d.y};
  struct vec r3 = {m.a.z, m.b.z, m.c.z, m.d.z};
  struct vec r4 = {m.a.a, m.b.a, m.c.a, m.d.a};
  set(&(p.a), dot(r1, n.a), dot(r2, n.a), dot(r3, n.a), dot(r4, n.a));
  set(&(p.b), dot(r1, n.b), dot(r2, n.b), dot(r3, n.b), dot(r4, n.b));
  set(&(p.c), dot(r1, n.c), dot(r2, n.c), dot(r3, n.c), dot(r4, n.c));
  set(&(p.d), dot(r1, n.d), dot(r2, n.d), dot(r3, n.d), dot(r4, n.d));
  return p;
}

// multiply a 4x4 matrix with a 4x1 vector to get another 4x1 vector
struct vec trans(const struct mat m, const struct vec v) {
  struct vec r1 = {m.a.x, m.b.x, m.c.x, m.d.x};
  struct vec r2 = {m.a.y, m.b.y, m.c.y, m.d.y};
  struct vec r3 = {m.a.z, m.b.z, m.c.z, m.d.z};
  struct vec r4 = {m.a.a, m.b.a, m.c.a, m.d.a};
  struct vec w = {dot(r1, v), dot(r2, v), dot(r3, v), dot(r4, v)};
  return w;
}

// get transpose of a matrix
struct mat transpose(const struct mat m) {
  struct vec r1 = {m.a.x, m.b.x, m.c.x, m.d.x};
  struct vec r2 = {m.a.y, m.b.y, m.c.y, m.d.y};
  struct vec r3 = {m.a.z, m.b.z, m.c.z, m.d.z};
  struct vec r4 = {m.a.a, m.b.a, m.c.a, m.d.a};
  // rows becom columns
  struct mat p = {r1, r2, r3, r4};
  return p;
}

// construct view matrix (like glm lookAt())
struct mat getViewMat(const struct vec eyePos, const struct vec lookPos,
    const struct vec upDir) {
  // n points behind camera
  struct vec n = normalize(vSub(eyePos, lookPos));
  // v points sideways
  struct vec v = normalize(cross(upDir, n));
  // u points up
  struct vec u = normalize(cross(n, v));

  struct vec temp = {0, 0, 0, 1};
  struct mat V = {v, u, n, temp};
  V = transpose(V);
  struct vec negRTd = {dot(scale(v, -1), eyePos), dot(scale(u, -1), eyePos),
      dot(scale(n, -1), eyePos), 1};
  V.d = negRTd;
  return V;
}