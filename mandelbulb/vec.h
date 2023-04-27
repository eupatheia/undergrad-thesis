#ifndef VEC_H_
#define VEC_H_

//////////////////////////////////////////////
// *** NOTE: ONLY dot() USES COMPONENT a ! ***
//////////////////////////////////////////////

// representation of a 3D vector/point (homogenous coordinates)
struct vec {
  float x;
  float y;
  float z;
  float a;  // mostly ignored for calculations
};

void set(struct vec * v, float x, float y, float z, float a);
void printVec(const struct vec v);

// vector arithmetic functions
struct vec vAdd(const struct vec u, const struct vec v);
struct vec vSub(const struct vec u, const struct vec v);
struct vec vMul(const struct vec u, const struct vec v);
struct vec vPow(const struct vec v, float exp);

float dot(const struct vec u, const struct vec v);
float distance(const struct vec u, const struct vec v);
float length(const struct vec v);
struct vec normalize(const struct vec v);
struct vec scale(const struct vec v, const float factor);
struct vec offset(const struct vec v, const float offset);
struct vec cross(const struct vec u, const struct vec v);
struct vec reflect(const struct vec v, const struct vec n);
struct vec refract(struct vec v, struct vec n, float refractiveIndexRatio);

#endif