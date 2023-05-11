#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "vec.h"

// clamp num to range [start, end]
float clampFloat(float num, float start, float end) {
  if (num > end) {
    return end;
  } else if (num < start) {
    return start;
  } else {
    return num;
  }
}

// clamp num to range [start, end]
int clampInt(int num, int start, int end) {
  if (num > end) {
    return end;
  } else if (num < start) {
    return start;
  } else {
    return num;
  }
}

// get density at localPos WRT center
float density(struct vec localPos) {
  // struct vec r = normalize(localPos);  // unit vector length 1
  // float rr = dot(r, r);
  // float pp = dot(localPos, localPos);
  // if (rr >= pp) {
  //   return 0.5;  // inside puff
  // }
  // float tmp = sqrt(rr) - (pp / rr);
  // return tmp > 0.000001 ? tmp : 0.0;
  float r = length(localPos);
  float coreRadius = 0.25;
  float coreDensity = 0.5;
  return clampFloat(exp(0.5 * -r), 0, 1);
}

// get color of puff at localPos WRT center
// adapted from https://www.shadertoy.com/view/wt33Wl
struct vec puffcolor(struct vec localPos) {
  struct vec inner = {0.8, 0.3, 0.0, 1.0};  
  struct vec outer = {0, 0.6, 0.47, 1.0};
  // use fraction of inner color based on distance to center
  float tmp = clampFloat(1.0 - length(localPos), 0, 1);
  // blend inner and outer color
  return vAdd(scale(inner, tmp), scale(outer, (1.0 - tmp)));
}

// return the face of the cube that contains pos:
// 0 = not on any face
// 1 = left = xmin
// 2 = right = xmax
// 3 = bottom = ymin
// 4 = top = ymax
// 5 = near = zmin
// 6 = far = zmax
int findFace(struct vec pos, float left, float right, float bottom, float top,
    float near, float far) {
  if (fabs(pos.x - left) < 0.1) {
    return 1;
  } else if (fabs(pos.x - right) < 0.1) {
    return 2;
  } else if (fabs(pos.y - bottom) < 0.1) {
    return 3;
  } else if (fabs(pos.y - top) < 0.1) {
    return 4;
  } else if (fabs(pos.z - near) < 0.1) {
    return 5;
  } else if (fabs(pos.z - far) < 0.1) {
    return 6;
  } else {
    return 0;  // not on any face
  }
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
  struct vec w = origPos;
  float x = w.x; float x2 = x*x; float x4 = x2*x2;
  float y = w.y; float y2 = y*y; float y4 = y2*y2;
  float z = w.z; float z2 = z*z; float z4 = z2*z2;

  float k3 = x2 + z2;
  float k2 = 1.0f / sqrt( k3*k3*k3*k3*k3*k3*k3 );
  float k1 = x4 + y4 + z4 - 6.0*y2*z2 - 6.0*x2*y2 + 2.0*z2*x2;
  float k4 = x2 - y2 + z2;

  w.x = origPos.x +  64.0*x*y*z*(x2-z2)*k4*(x4-6.0*x2*z2+z4)*k1*k2;
  w.y = origPos.y + -16.0*y2*k3*k4*k4 + k1*k1;
  w.z = origPos.z +  -8.0*y*k4*(x4*x4 - 28.0*x4*x2*z2 + 70.0*x4*z4 - 28.0*x2*z2*z4 + z4*z4)*k1*k2;
  return w;
}

// calculate SDF (signed distance function) from a position pos
// to the mandelbulb object, returning distance
// adapted from https://iquilezles.org/articles/mandelbulb/:
float sdScene(struct vec pos, int maxIterations) {
  struct vec w = pos;
  float dw = 1.0;  // gradient of potential
  float m = length(w);  // modulus of a point (i.e. length, |w|)
  int iter = 0;
  while (iter < maxIterations && m < 16) {
    // compute gradient dw_{n+1} = 8 * |w_n|^7 * dw_n
    dw = 8.0 * pow(m, 7) * dw + 1.0;
    // get next iteration
    w = optf(pos, w);
    m = length(w);
    iter++;
  }
  // calc distance d = |w| * log|w| / |dw|;
  // note this is modulus, not abs val;
  // multiply 0.25 to offset errors in sdf
  // (see https://iquilezles.org/articles/distancefractals/)
  return (m * log(m) / dw) * 0.25;
}


// calculate SDF of a position pos and direction ray from an axis-aligned
// bounding box specified by [xmin, xmax], [ymin, ymax], [zmin, zmax]
float sdBox(struct vec pos, struct vec ray, float xmin, float xmax, float ymin,
    float ymax, float zmin, float zmax) {
  float tmin = -FLT_MAX;
  float tmax = FLT_MAX;
  float t1, t2;
  if (ray.x != 0) {
    t1 = (xmin - pos.x) / ray.x;
    t2 = (xmax - pos.x) / ray.x;
    tmin = fmax(tmin, fmin(t1, t2));
    tmax = fmin(tmax, fmax(t1, t2));
  }
  if (ray.y != 0) {
    t1 = (ymin - pos.y) / ray.y;
    t2 = (ymax - pos.y) / ray.y;
    tmin = fmax(tmin, fmin(t1, t2));
    tmax = fmin(tmax, fmax(t1, t2));
  }
  if (ray.z != 0) {
    t1 = (zmin - pos.z) / ray.z;
    t2 = (zmax - pos.z) / ray.z;
    tmin = fmax(tmin, fmin(t1, t2));
    tmax = fmin(tmax, fmax(t1, t2));
  }
  if (tmin <= tmax) {
    if (tmin > 0) {
      return tmin;
    } else {
      // ray starts inside box, tmin is behind
      return tmax;
    }
  } else {
    return -1;  // no intersection
  }
}

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

// calculate the normal at some position pos on a surface
// by calculating the gradient in every direction (how does a change
// by epsilon along one axis change the SDF?)
struct vec calcNormal(struct vec pos, float epsilon, int maxIterations) {
  struct vec x1 = {pos.x - epsilon, pos.y, pos.z, 1.0};
  struct vec x2 = {pos.x + epsilon, pos.y, pos.z, 1.0};
  struct vec y1 = {pos.x, pos.y - epsilon, pos.z, 1.0};
  struct vec y2 = {pos.x, pos.y + epsilon, pos.z, 1.0};
  struct vec z1 = {pos.x, pos.y, pos.z - epsilon, 1.0};
  struct vec z2 = {pos.x, pos.y, pos.z + epsilon, 1.0};
  float deltaX = sdScene(x2, maxIterations) - sdScene(x1, maxIterations);
  float deltaY = sdScene(y2, maxIterations) - sdScene(y1, maxIterations);
  float deltaZ = sdScene(z2, maxIterations) - sdScene(z1, maxIterations);
  struct vec norm = {deltaX, deltaY, deltaZ, 0.0};
  return normalize(norm);  // return unit vector
}

// distance-aided raymarch:
// given a normalized ray with tail at pos,
// performs raymarch and returns 1 if hit, else 0,
// and if hit, returns position of surface in hitPos
int DAraymarch(struct vec pos, struct vec ray, struct vec * hitPos,
    float * shadowColor, int maxSteps, float hitRange, int maxIterations,
    float softness) {
  // compute distance to bounding sphere to only search w/in bounds
  struct vec center = {0, 0, 0, 0};
  struct vec sphereDistance = sphereIntersect(pos, ray, center, 1.25f);
  if (sphereDistance.y < 0) {
    // no far intersection, never hit bounding sphere, never hit mandelbulb
    return 0;
  }
  // start at pos or near intersection (but not behind pos)
  float totalDistance = fmax(0, sphereDistance.x);
  int steps = 0;
  *shadowColor = 1;  // start with full brightness
  // stop when too many steps or went past far intersection
  while (steps < maxSteps && totalDistance <= sphereDistance.y) {
    struct vec currPos = vAdd(pos, scale(ray, totalDistance));
    float minSDF = fabsf(sdScene(currPos, maxIterations));
    if (minSDF < hitRange) {
      // close enough, call a hit
      hitPos->x = currPos.x;
      hitPos->y = currPos.y;
      hitPos->z = currPos.z;
      hitPos->a = currPos.a;
      *shadowColor = 0.2;
      return 1;
    }
    totalDistance += minSDF;
    steps++;
    // clamp between 0.1 and 1
    *shadowColor = fmax(fmin(*shadowColor, softness * minSDF / totalDistance), 0.2);
  }
  // never hit mandelbulb surface in scene
  return 0;
}

// REGULAR equal-step raymarch:
// given a normalized ray with tail at pos,
// performs raymarch and returns 1 if hit, else 0,
// and if hit, returns position of surface in hitPos
int STEPraymarch(struct vec pos, struct vec ray, struct vec * hitPos,
    float maxMarchDistance, float stepSize, int maxIterations) {
  float d = 0;
  for (float i = 0; i < maxMarchDistance; i += stepSize) {
    struct vec currPos = vAdd(pos, scale(ray, i));
    float d = sdScene(currPos, maxIterations);
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

// return far side bulb intersection along ray
struct vec farBulbIntersect(struct vec origPos, struct vec ray,
    struct vec hitPos, int maxIterations, int maxSteps, float hitRange,
    float softness) {
  ray = normalize(ray);
  struct vec center = {0, 0, 0, 1};
  struct vec intersections = sphereIntersect(origPos, ray, center, 1.25);
  struct vec oppositePos = vAdd(origPos, scale(ray, intersections.y));
  float shadowColor;
  int hit = DAraymarch(oppositePos, scale(ray, -1), &hitPos, &shadowColor,
      maxSteps, hitRange, maxIterations, softness);
  return hitPos;
}