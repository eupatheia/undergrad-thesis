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

struct vec vPow(const struct vec v, float exp) {
  struct vec w;
  w.x = pow(v.x, exp);
  w.y = pow(v.y, exp);
  w.z = pow(v.z, exp);
  w.a = 0.0f;
  return w;
}

float dot(const struct vec u, const struct vec v) {
  return (u.x * v.x) + (u.y * v.y) + (u.z * v.z);
}

float distance(const struct vec u, const struct vec v) {
  return sqrt(pow(u.x - v.x, 2) + pow(u.y - v.y, 2) + pow(u.z - v.z, 2));
}

float length(const struct vec v) {
  struct vec zero = {0, 0, 0, 0};
  return distance(v, zero);
}

// if treating v as vector, set v.a = 0.0 before passing to this function,
// i.e. should not be called if v is treated as a point (v.a = 1.0)
struct vec normalize(const struct vec v) {
  float len = sqrt(pow(v.x, 2) + pow(v.y, 2) + pow(v.z, 2));
  return scale(v, 1.0f / len);
}

// multiply a vector by a scaling factor
struct vec scale(const struct vec v, const float factor) {
  struct vec w;
  w.x = v.x * factor;
  w.y = v.y * factor;
  w.z = v.z * factor;
  w.a = 0.0f;
  return w;
}

// add a number to vector, component-wise
struct vec offset(const struct vec v, const float offset) {
  struct vec w;
  w.x = v.x + offset;
  w.y = v.y + offset;
  w.z = v.z + offset;
  w.a = 0.0f;
  return w;
}