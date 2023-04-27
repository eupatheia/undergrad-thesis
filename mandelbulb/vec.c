#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "vec.h"

// reset values inside the struct
void set(struct vec * v, float x, float y, float z, float a) {
  v->x = x;
  v->y = y;
  v->z = z;
  v->a = a;
}

void printVec(const struct vec v) {
  printf("(%.3f %.3f %.3f %.3f)\n", v.x, v.y, v.z, v.a);
}

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
  return (u.x * v.x) + (u.y * v.y) + (u.z * v.z) + (u.a * v.a);
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

// cross product of two vectors
struct vec cross(const struct vec u, const struct vec v) {
  struct vec w;
  w.x = (u.y * v.z) - (u.z * v.y);
  w.y = (u.z * v.x) - (u.x * v.z);
  w.z = (u.x * v.y) - (u.y * v.x);
  w.a = 1.0f;
  return w;
}

// return the ray resulting from reflecting v about normal n
struct vec reflect(const struct vec v, const struct vec n) {
  return vAdd(scale(n, -2 * dot(v, n)), v);
}

// return the ray resulting from refracting v off a surface with normal n,
// given a ratio of refractive indices for each color channel
struct vec refract(struct vec v, struct vec n,
    float refractiveIndexRatio) {
  v = normalize(v);
  n = normalize(n);
  float cosTheta = dot(scale(v, -1), n) / (length(v) * length(n));
  struct vec perpendicular = scale(vAdd(v, scale(n, cosTheta)),
      refractiveIndexRatio);
  float modPerpSq = pow(length(perpendicular), 2);
  if (modPerpSq > 1) {
    // total internal reflection, no refraction
    return reflect(v, n);
  } else {
    // struct vec parallel = scale(n, -sqrt(1 - (refractiveIndexRatio *
    //     refractiveIndexRatio * (1 - (cosTheta * cosTheta)))));
    struct vec parallel = scale(n, -sqrt(1 - modPerpSq));
    return vAdd(perpendicular, parallel);
  }
}