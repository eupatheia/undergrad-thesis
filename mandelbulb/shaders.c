#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "read_write.h"
#include "vec.h"
#include "mat.h"
#include "shaders.h"
#include "ray_functions.h"

// calculate phong shading colors
struct lighting phongShading(struct vec pos, struct vec norm, struct vec La,
    struct vec Ld, struct vec Ls, struct vec Ka, struct vec Kd,
    struct vec Ks, float shininess, struct vec lightPos,
    struct vec camPos) {
  // normalized direction to light source
  struct vec s = normalize(vSub(lightPos, pos));
  // normalized direction to viewer
  struct vec v = {camPos.x - pos.x, camPos.y - pos.y,
      camPos.z - pos.z, 0.0f};
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

// return background color based on ratio of y position to image height,
// yratio in [0, 1]
struct ppm_pixel backgroundColor(float yratio) {
  struct ppm_pixel color;
  // clamp to [0, 1]
  yratio = fmax(fmin(yratio, 1), 0);
  color.red = yratio * 255;
  color.green = yratio * 255;
  color.blue = yratio * 255;
  return color;
}

// return color from point with ypos on face of a cube,
// assuming top and bottom have specific colors and sides interpolate
struct ppm_pixel twoPlaneCubemap(int face, float ypos) {
  struct ppm_pixel color;
  if (face == 3) {
    // color of bottom face
    color.red = 200;
    color.green = 200;
    color.blue = 200;
    return color;
  } else if (face == 4) {
    // color of top face
    color.red = 50;
    color.green = 50;
    color.blue = 50;
    return color;
  } else { // face == 1 || face == 2 || face == 5 || face == 6
    float t = (ypos + 10) / 20.0;
    // color = bottomColor * (1-t) + topColor * t
    color.red = 200 * (1 - t) + 50 * t;
    color.green = 200 * (1 - t) + 50 * t;
    color.blue = 200 * (1 - t) + 50 * t;
    return color;
    // return backgroundColor(t);
  }
}

// color by normals
struct ppm_pixel normalShader(float yratio, struct vec norm) {
  struct ppm_pixel color;
  color.red = (norm.x + 1) * 0.5 * 255;
  color.green = (norm.y + 1) * 0.5 * 255;
  color.blue = (norm.z + 1) * 0.5 * 255;
  return color;
}

// phong shader with soft shadows
struct ppm_pixel phongShadowShader(struct vec pos, struct vec hitPos,
    struct vec norm, struct vec lightPos, float hitRange, int maxSteps,
    int maxIterations, float softness) {
  // calculate phong shading color first
  struct vec lightColor = {1.0, 1.0, 1.0, 0.0};
  struct vec Ka = {0.24725f, 0.2245f, 0.0645f, 0.0f};
  struct vec Kd = {0.34615f, 0.3143f, 0.0903f, 0.0f};
  struct vec Ks = {0.797357f, 0.72399f, 0.20801f, 0.0f};
  // struct vec Ka = {0, 0, 0.2f, 0};
  // struct vec Kd = {0, 0.3, 0.7, 0};
  // struct vec Ks = {1, 1, 1, 0};
  struct lighting shades = phongShading(hitPos, norm, lightColor, lightColor,
      lightColor, Ka, Kd, Ks, 20.0f, lightPos, pos);

  // normalized direction to point light source
  // struct vec lightDir = normalize(vSub(lightPos, hitPos));
  // normalized direction to directional light source
  struct vec lightDir = {-5, 5, 5, 0};
  lightDir = normalize(lightDir);
  struct vec shadowPos = {0, 0, 0, 0};
  // start slightly away from surface to avoid hitting same point again
  struct vec hitPosPlus = vAdd(hitPos, scale(lightDir, hitRange * 2));
  // then calculate shadow ray
  float shadowColor;
  int hasShadow = DAraymarch(hitPosPlus, lightDir, &shadowPos, &shadowColor,
      maxSteps, hitRange, maxIterations, softness);
  struct vec finalColor = scale(vAdd(vAdd(shades.ambient, shades.diffuse),
      shades.specular), shadowColor);
  struct ppm_pixel color;
  color.red = fmin(finalColor.x * 255, 255);
  color.green = fmin(finalColor.y * 255, 255);
  color.blue = fmin(finalColor.z * 255, 255);
  return color;
}

// simulate reflection color of ray about norm at hitPos
struct ppm_pixel reflectShader(struct vec ray, struct vec norm,
    struct vec hitPos) {
  // reflect viewing ray about the surface normal
  struct vec r = normalize(reflect(ray, norm));
  float distanceToBox = sdBox(hitPos, r, -10, 10, -10, 10, -10, 10);
  struct vec cubePoint = vAdd(hitPos, scale(r, distanceToBox));
  int face = findFace(cubePoint, -10, 10, -10, 10, -10, 10);
  if (face != 0) {
    return twoPlaneCubemap(face, cubePoint.y);
  } else { // face == 0
    // must intersect a face, error
    printf("Error, ray did not intersect the skybox: ");
    printVec(cubePoint);
    printVec(hitPos);
    printVec(ray);
    printVec(norm);
    printVec(r);
    printf("%.3f\n", distanceToBox);
    exit(0);
  }
}

// simulate refraction color of ray about norm at hitPos,
// through materials with refractive index ratio
struct ppm_pixel refractShader(struct vec ray, struct vec norm,
    struct vec hitPos, float refractiveIndexRatio) {
  // refract viewing ray at the surface normal
  struct vec r = refract(ray, norm, refractiveIndexRatio);
  float distanceToBox = sdBox(hitPos, r, -10, 10, -10, 10, -10, 10);
  struct vec cubePoint = vAdd(hitPos, scale(r, distanceToBox));
  int face = findFace(cubePoint, -10, 10, -10, 10, -10, 10);
  if (face != 0) {
    return twoPlaneCubemap(face, cubePoint.y);
  } else { // face == 0
    // must intersect a face, error
    printf("Error, ray did not intersect the skybox: ");
    printVec(cubePoint);
    printVec(hitPos);
    printVec(ray);
    printVec(norm);
    printVec(r);
    printf("%.3f\n", distanceToBox);
    exit(0);
  }
}

// dielectric shader (combo of reflection and refraction colors using 
// fresnel coefficient for mixing)
struct ppm_pixel dielectricShader(struct vec ray, struct vec norm,
    struct vec hitPos, float refractiveIndexRatio) {
  ray = normalize(ray);
  norm = normalize(norm);
  // clamp fresnel between [0, 1]
  float fresnel = fmax(fmin(0.1 + (2 * pow(1 + dot(ray, norm), 2)), 1), 0);
  struct ppm_pixel reflectColor = reflectShader(ray, norm, hitPos);
  struct ppm_pixel refractColor = refractShader(ray, norm, hitPos,
      refractiveIndexRatio);
  struct ppm_pixel color;
  color.red = ((1 - fresnel) * reflectColor.red) + (fresnel * refractColor.red);
  color.green = ((1 - fresnel) * reflectColor.green) + (fresnel * refractColor.green);
  color.blue = ((1 - fresnel) * reflectColor.blue) + (fresnel * refractColor.blue);
  return color;
}

// chromatic dispersion shader (different refractiveIndexRatio for each RGB)
struct ppm_pixel chromaticDispersionShader(struct vec ray, struct vec norm,
    struct vec hitPos, struct vec refractiveIndexRatios) {
  struct ppm_pixel color;
  color.red = refractShader(ray, norm, hitPos, refractiveIndexRatios.x).red;
  color.green = refractShader(ray, norm, hitPos, refractiveIndexRatios.y).green;
  color.blue = refractShader(ray, norm, hitPos, refractiveIndexRatios.z).blue;
  return color;
}

// simulate volume by step-raymarching through;
// adapted from https://www.shadertoy.com/view/wt33Wl
struct ppm_pixel volumetricShader(struct vec origPos, struct vec ray,
    struct vec hitPos, struct vec norm, float yratio, int maxIterations,
    float refractiveIndexRatio, int maxSteps, float hitRange,
    float softness) {
  ray = normalize(ray);
  struct vec center = {0, 0, 0, 0};  // center of object
  float transmitivity = 1.0;
  float kappa = 0.2;  // hyperparam
  struct vec c = {0.0, 0.0, 0.0, 1.0};  // accumulated color
  float stepSize = 0.01;
  // calculate far intersection with bulb along ray
  struct vec farHitPos = farBulbIntersect(origPos, ray, hitPos,
      maxIterations, maxSteps, hitRange, softness);
  float range = length(vSub(farHitPos, hitPos));
  // printf("%.3f\n", range);
  for (float t = 0; t < range; t += stepSize) {
    // step along ray, currPos is already local WRT center at origin
    struct vec currPos = vAdd(hitPos, scale(ray, t));
    // calculate transmitivity based on current density
    float d = density(currPos);
    float transmitivityDelta = exp(-kappa * stepSize * d);
    transmitivity *= transmitivityDelta;
    struct vec pcolor = puffcolor(currPos);
    // accumulate color from current position
    c = vAdd(c, scale(pcolor,
        ((1.0 - transmitivityDelta) / kappa) * transmitivity));
  }
  struct ppm_pixel background = backgroundColor(yratio);
  struct ppm_pixel color;
  // get final color and clamp to range
  color.red = clampInt(((background.red / 255.0 * transmitivity) + c.x)
      * 255, 0, 255);
  color.green = clampInt(((background.green / 255.0 * transmitivity) + c.y)
      * 255, 0, 255);
  color.blue = clampInt(((background.blue / 255.0 * transmitivity) + c.z)
      * 255, 0, 255);
  return color;
}

// blend dielectric and volumetric shading
struct ppm_pixel blendShader(struct vec origPos, struct vec ray,
    struct vec hitPos, struct vec norm, float yratio, int maxIterations,
    float refractiveIndexRatio, int maxSteps, float hitRange,
    float softness) {
  ray = normalize(ray);
  struct vec center = {0, 0, 0, 0};  // center of object
  float transmitivity = 1.0;
  float kappa = 0.2;  // hyperparam
  struct vec c = {0.0, 0.0, 0.0, 1.0};  // accumulated color
  float stepSize = 0.01;
  // calculate far intersection with bulb along ray
  struct vec farHitPos = farBulbIntersect(origPos, ray, hitPos,
      maxIterations, maxSteps, hitRange, softness);
  float range = length(vSub(farHitPos, hitPos));
  // printf("%.3f\n", range);
  for (float t = 0; t < range; t += stepSize) {
    // step along ray, currPos is already local WRT center at origin
    struct vec currPos = vAdd(hitPos, scale(ray, t));
    // calculate transmitivity based on current density
    float d = density(currPos);
    float transmitivityDelta = exp(-kappa * stepSize * d);
    transmitivity *= transmitivityDelta;
    struct vec pcolor = puffcolor(currPos);
    // accumulate color from current position
    c = vAdd(c, scale(pcolor,
        ((1.0 - transmitivityDelta) / kappa) * transmitivity));
  }
  struct ppm_pixel background = dielectricShader(ray, norm, hitPos, 0.9);
  struct ppm_pixel color;
  // get final color and clamp to range
  color.red = clampInt(((background.red / 255.0 * transmitivity) + c.x)
      * 255, 0, 255);
  color.green = clampInt(((background.green / 255.0 * transmitivity) + c.y)
      * 255, 0, 255);
  color.blue = clampInt(((background.blue / 255.0 * transmitivity) + c.z)
      * 255, 0, 255);
  return color;
}