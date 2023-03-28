#ifndef VEC_H_
#define VEC_H_

// representation of a 3D vector/point (homogenous coordinates)
struct vec {
  float x;
  float y;
  float z;
  float a;  // mostly ignored for calculations
};

// vector arithmetic functions
struct vec vAdd(const struct vec u, const struct vec v);
struct vec vSub(const struct vec u, const struct vec v);
struct vec vMul(const struct vec u, const struct vec v);
struct vec vDiv(const struct vec u, const struct vec v);

float dot(const struct vec u, const struct vec v);
float distance(const struct vec u, const struct vec v);
struct vec normalize(const struct vec v);

#endif