#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "vec.h"

struct vec vAdd(const struct vec u, const struct vec v) {
  struct vec w;
  w.x = u.x + v.x;
  w.y = u.y + v.y;
  w.z = u.z + v.z;
  w.a = 0.0f;
  return w;
}

struct vec vSub(const struct vec u, const struct vec v) {
  struct vec w;
  w.x = u.x - v.x;
  w.y = u.y - v.y;
  w.z = u.z - v.z;
  w.a = 0.0f;
  return w;
}

struct vec vMul(const struct vec u, const struct vec v) {
  struct vec w;
  w.x = u.x * v.x;
  w.y = u.y * v.y;
  w.z = u.z * v.z;
  w.a = 0.0f;
  return w;
}

struct vec vDiv(const struct vec u, const struct vec v) {
  struct vec w;
  w.x = u.x / v.x;
  w.y = u.y / v.y;
  w.z = u.z / v.z;
  w.a = 0.0f;
  return w;
}

float dot(const struct vec u, const struct vec v) {
  return (u.x * v.x) + (u.y * v.y) + (u.z * v.z);
}

float distance(const struct vec u, const struct vec v) {
  return sqrt(pow(u.x - v.x, 2) + pow(u.y - v.y, 2) + pow(u.z - v.z, 2));
}

// if treating v as vector, set v.a = 0.0 before passing to this function,
// i.e. should not be called if v is treated as a point (v.a = 1.0)
struct vec normalize(const struct vec v) {
  float len = sqrt(pow(v.x, 2) + pow(v.y, 2) + pow(v.z, 2));
  struct vec w = {len, len, len, len};
  return vDiv(v, w);
}