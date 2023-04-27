// functions that calulate color of object

#ifndef SHADERS_H_
#define SHADERS_H_

#include "read_write.h"
#include "vec.h"

struct lighting {
  struct vec ambient;
  struct vec diffuse;
  struct vec specular;
};

struct lighting phongShading(struct vec pos, struct vec norm, struct vec La,
    struct vec Ld, struct vec Ls, struct vec Ka, struct vec Kd,
    struct vec Ks, float shininess, struct vec lightPos,
    struct vec camPos);
struct ppm_pixel backgroundColor(float yratio);
struct ppm_pixel twoPlaneCubemap(int face, float ypos);

struct ppm_pixel normalShader(float yratio, struct vec norm);
struct ppm_pixel phongShadowShader(struct vec pos, struct vec hitPos,
    struct vec norm, struct vec lightPos, float hitRange, int maxSteps,
    int maxIterations, float softness);
struct ppm_pixel reflectShader(struct vec ray, struct vec norm,
    struct vec hitPos);
struct ppm_pixel refractShader(struct vec ray, struct vec norm,
    struct vec hitPos, float refractiveIndexRatio);
struct ppm_pixel dielectricShader(struct vec ray, struct vec norm,
    struct vec hitPos, float refractiveIndexRatio);
struct ppm_pixel chromaticDispersionShader(struct vec ray, struct vec norm,
    struct vec hitPos, struct vec refractiveIndexRatios);
struct ppm_pixel volumetricShader(struct vec origPos, struct vec ray,
    struct vec hitPos, struct vec norm, float yratio, int maxIterations,
    float refractiveIndexRatio, int maxSteps, float hitRange,
    float softness);
#endif