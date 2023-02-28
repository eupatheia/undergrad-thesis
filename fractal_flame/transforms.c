#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include "transforms.h"
#include "read_write.h"

// perform basic affine transformation on point p,
// i.e. F(x, y) = (ax + by + c, dx + ey + f);
extern struct point affine(struct point p, float a, float b, float c, float d,
    float e, float f) {
  struct point next;
  next.x = (a * p.x) + (b * p.y) + c;
  next.y = (d * p.x) + (e * p.y) + f;
  return next;
}

// rotate point p by angle (in radians)
extern struct point rotate(struct point p, float angle) {
  struct point next;
  next.x = (p.x * cos(angle)) - (p.y * sin(angle));
  next.y = (p.y * cos(angle)) + (p.x * sin(angle));
  return next;
}

// generate random float between start and end, inclusive
extern float randomParam(float start, float end) {
  return (((float) rand() / RAND_MAX) * (end - start)) + start;
}

// =========================================================
// =========================================================

extern struct point linear(struct point p, float a, float b, float c, float d,
    float e, float f) {
  p = affine(p, a, b, c, d, e, f);
  return p;  // no change
}
extern struct point sinusoidal(struct point p, float a, float b, float c,
    float d, float e, float f) {
  p = affine(p, a, b, c, d, e, f);
  p.x = sin(p.x);
  p.y = sin(p.y);
  return p;
}

extern struct point spherical(struct point p, float a, float b, float c,
    float d, float e, float f) {
  p = affine(p, a, b, c, d, e, f);
  float r = sqrt(pow(p.x, 2) + pow(p.y, 2));
  struct point next;
  next.x = p.x / pow(r, 2);
  next.y = p.y / pow(r, 2);
  return next;
}
extern struct point swirl(struct point p, float a, float b, float c,
    float d, float e, float f) {
  p = affine(p, a, b, c, d, e, f);
  float r = sqrt(pow(p.x, 2) + pow(p.y, 2));
  struct point next;
  next.x = (p.x * sin(pow(r, 2))) - (p.y * cos(pow(r, 2)));
  next.y = (p.x * cos(pow(r, 2))) + (p.y * sin(pow(r, 2)));
  return next;
}

extern struct point horseshoe(struct point p, float a, float b, float c, float d,
    float e, float f) {
  p = affine(p, a, b, c, d, e, f);
  float r = sqrt(pow(p.x, 2) + pow(p.y, 2));
  struct point next;
  next.x = ((p.x - p.y) * (p.x + p.y)) / r;
  next.y = (2 * p.x * p.y) / r;
  return next;
}

extern struct point polar(struct point p, float a, float b, float c,
    float d, float e, float f) {
  p = affine(p, a, b, c, d, e, f);
  float r = sqrt(pow(p.x, 2) + pow(p.y, 2));
  float theta = atan((float) p.x / p.y);
  struct point next;
  next.x = theta / M_PI;
  next.y = r - 1;
  return next;
}

extern struct point handkerchief(struct point p, float a, float b, float c,
    float d, float e, float f) {
  p = affine(p, a, b, c, d, e, f);
  float r = sqrt(pow(p.x, 2) + pow(p.y, 2));
  float theta = atan((float) p.x / p.y);
  struct point next;
  next.x = r * sin(theta + r);
  next.y = r * cos(theta - r);
  return next;
}

extern struct point heart(struct point p, float a, float b, float c,
    float d, float e, float f) {
  p = affine(p, a, b, c, d, e, f);
  float r = sqrt(pow(p.x, 2) + pow(p.y, 2));
  float theta = atan((float) p.x / p.y);
  struct point next;
  next.x = r * sin(theta * r);
  next.y = r * -cos(theta * r);
  return next;
}

extern struct point disk(struct point p, float a, float b, float c,
    float d, float e, float f) {
  p = affine(p, a, b, c, d, e, f);
  float r = sqrt(pow(p.x, 2) + pow(p.y, 2));
  float theta = atan((float) p.x / p.y);
  struct point next;
  next.x = (theta / M_PI) * sin(M_PI * r);
  next.y = (theta / M_PI) * cos(M_PI * r);
  return next;
}

extern struct point spiral(struct point p, float a, float b, float c,
    float d, float e, float f) {
  p = affine(p, a, b, c, d, e, f);
  float r = sqrt(pow(p.x, 2) + pow(p.y, 2));
  float theta = atan((float) p.x / p.y);
  struct point next;
  next.x = (cos(theta) + sin(r)) / r;
  next.y = (sin(theta) - cos(r)) / r;
  return next;
}

extern struct point hyperbolic(struct point p, float a, float b, float c,
    float d, float e, float f) {
  p = affine(p, a, b, c, d, e, f);
  float r = sqrt(pow(p.x, 2) + pow(p.y, 2));
  float theta = atan((float) p.x / p.y);
  struct point next;
  next.x = sin(theta) / r;
  next.y = r * cos(theta);
  return next;
}

extern struct point diamond(struct point p, float a, float b, float c,
    float d, float e, float f) {
  p = affine(p, a, b, c, d, e, f);
  float r = sqrt(pow(p.x, 2) + pow(p.y, 2));
  float theta = atan((float) p.x / p.y);
  struct point next;
  next.x = sin(theta) * cos(r);
  next.y = cos(theta) * sin(r);
  return next;
}

extern struct point ex(struct point p, float a, float b, float c,
    float d, float e, float f) {
  p = affine(p, a, b, c, d, e, f);
  float r = sqrt(pow(p.x, 2) + pow(p.y, 2));
  float theta = atan((float) p.x / p.y);
  float p0 = sin(theta + r);
  float p1 = cos(theta - r);
  struct point next;
  next.x = r * (pow(p0, 3) + pow(p1, 3));
  next.y = r * (pow(p0, 3) - pow(p1, 3));
  return next;
}

extern struct point julia(struct point p, float a, float b, float c,
    float d, float e, float f) {
  p = affine(p, a, b, c, d, e, f);
  float r = sqrt(pow(p.x, 2) + pow(p.y, 2));
  float theta = atan((float) p.x / p.y);
  float omega;
  float k = rand() % 2;
  if (k == 0) {
    omega = 0;
  } else {
    omega = M_PI;
  }
  struct point next;
  next.x = sqrt(r) * cos((theta / 2) + omega);
  next.y = sqrt(r) * sin((theta / 2) + omega);
  return next;
}

extern struct point fisheye(struct point p, float a, float b, float c,
    float d, float e, float f) {
  p = affine(p, a, b, c, d, e, f);
  float r = sqrt(pow(p.x, 2) + pow(p.y, 2));
  struct point next;
  next.x = (2 / (r + 1)) * p.y;
  next.y = (2 / (r + 1)) * p.x;
  return next;
}

extern struct point exponential(struct point p, float a, float b, float c,
    float d, float e, float f) {
  p = affine(p, a, b, c, d, e, f);
  float r = sqrt(pow(p.x, 2) + pow(p.y, 2));
  struct point next;
  next.x = exp(p.x - 1) * cos(M_PI * p.y);
  next.y = exp(p.x - 1) * sin(M_PI * p.y);
  return next;
}

// // given parameter val
// extern struct point rings2(struct point p, float val, float a, float b, float c,
//     float d, float e, float f) {
//   p = affine(p, a, b, c, d, e, f);
//   float r = sqrt(pow(p.x, 2) + pow(p.y, 2));
//   float theta = atan((float) p.x / p.y);
//   float param = pow(val, 2);
//   float t = r - (2 * param * floor((r + param) / (2 * param))) +
//       (r * (1 - param));
//   struct point next;
//   next.x = t * sin(theta);
//   next.y = t * cos(theta);
//   return next;
// }

extern struct point eyefish(struct point p, float a, float b, float c,
    float d, float e, float f) {
  p = affine(p, a, b, c, d, e, f);
  float r = sqrt(pow(p.x, 2) + pow(p.y, 2));
  struct point next;
  next.x = (2 / (r + 1)) * p.x;
  next.y = (2 / (r + 1)) * p.y;
  return next;
}

// =========================================================
// =========================================================

extern void systemSphericalSpherical(struct point * p, float * c) {
  // associated color for F_0
  float c0 = 1.0;
  // associated color for F_1
  float c1 = 0.0;
  int k = rand() % 2;  // both have probability 0.5, here
  if (k == 0) {
    *p = spherical(*p, 0.562482, -0.539599, -0.42992, 0.397861, 0.501088, -0.112404);
    *c = (*c + c0) / 2;
  } else {  // k == 1
    *p = spherical(*p, 0.830039, 0.16248, 0.91022, -0.496174, 0.75046, 0.288389);
    *c = (*c + c1) / 2;
  }
}

extern void systemLinearSwirlSpiral(struct point * p, float * c) {
  // associated colors per function
  float c0 = 1.0;
  float c1 = 0.0;
  float c2 = 0.5;
  int k = rand() % 3;
  if (k == 0) {
    *p = linear(*p, 0.98396, 0.359416, -0.85059, 0.298328, -0.583541, -0.378754);
    *c = (*c + c0) / 2;
  } else if (k == 1) {
    *p = spiral(*p, -0.900388, 0.397598, 0.465126, 0.293063, 0.0225629, -0.277212);
    *c = (*c + c1) / 2;
  } else {  // k == 2
    *p = swirl(*p, -0.329863, -0.369381, 0.977861, -0.0855261, -0.858379, 0.547595);
    *c = (*c + c2) / 2;
  }
}

extern void systemLinearSwirlSpiral_3SYM(struct point * p, float * c) {
  // associated colors per function
  float c0 = 1.0;
  float c1 = 0.0;
  float c2 = 0.5;
  int k = rand() % 9;
  if (k < 3) {
    // 120 degree rotation, no change in color to maintain symmetry
    *p = rotate(*p, (120 * M_PI) / 180);
  } else if (k < 6) {
    // 240 degree rotation, no change in color to maintain symmetry
    *p = rotate(*p, (240 * M_PI) / 180);
  } else if (k == 6) {
    *p = linear(*p, 0.98396, 0.359416, -0.85059, 0.298328, -0.583541, -0.378754);
    *c = (*c + c0) / 2;
  } else if (k == 7) {
    *p = spiral(*p, -0.900388, 0.397598, 0.465126, 0.293063, 0.0225629, -0.277212);
    *c = (*c + c1) / 2;
  } else {  // k == 8
    *p = swirl(*p, -0.329863, -0.369381, 0.977861, -0.0855261, -0.858379, 0.547595);
    *c = (*c + c2) / 2;
  }
}

extern void pickFrom2(struct point * p, float * col, transform_ptr * functions,
    float * weights, float * affineParams, float k, float booster) {
  // associated colors per function
  float col0 = 1.0;
  float col1 = 0.0;
  if (k < booster + weights[0]) {
    *p = functions[0](*p, affineParams[0], affineParams[1], affineParams[2],
        affineParams[3], affineParams[4], affineParams[5]);
    *col = (*col + col0) / 2;
  } else {
    *p = functions[1](*p, affineParams[6], affineParams[7], affineParams[8],
        affineParams[9], affineParams[10], affineParams[11]);
    *col = (*col + col1) / 2;
  }
}

extern void pickFrom3(struct point * p, float * col, transform_ptr * functions,
    float * weights, float * affineParams, float k, float booster) {
  // associated colors per function
  float col0 = 1.0;
  float col1 = 0.0;
  float col2 = 0.5;
  if (k < booster + weights[0]) {
    *p = functions[0](*p, affineParams[0], affineParams[1], affineParams[2],
        affineParams[3], affineParams[4], affineParams[5]);
    *col = (*col + col0) / 2;
  } else if (k < booster + weights[0] + weights[1]) {
    *p = functions[1](*p, affineParams[6], affineParams[7], affineParams[8],
        affineParams[9], affineParams[10], affineParams[11]);
    *col = (*col + col1) / 2;
  } else {
    *p = functions[2](*p, affineParams[12], affineParams[13], affineParams[14],
        affineParams[15], affineParams[16], affineParams[17]);
    *col = (*col + col2) / 2;
  }
}

extern void pickFrom4(struct point * p, float * col, transform_ptr * functions,
    float * weights, float * affineParams, float k, float booster) {
  // associated colors per function
  float col0 = 1.0;
  float col1 = 0.0;
  float col2 = 0.33;
  float col3 = 0.67;
  if (k < booster + weights[0]) {
    *p = functions[0](*p, affineParams[0], affineParams[1], affineParams[2],
        affineParams[3], affineParams[4], affineParams[5]);
    *col = (*col + col0) / 2;
  } else if (k < booster + weights[0] + weights[1]) {
    *p = functions[1](*p, affineParams[6], affineParams[7], affineParams[8],
        affineParams[9], affineParams[10], affineParams[11]);
    *col = (*col + col1) / 2;
  } else if (k < booster + weights[0] + weights[1] + weights[2]) {
    *p = functions[2](*p, affineParams[12], affineParams[13], affineParams[14],
        affineParams[15], affineParams[16], affineParams[17]);
    *col = (*col + col2) / 2;
  } else {
    *p = functions[3](*p, affineParams[18], affineParams[19], affineParams[20],
        affineParams[21], affineParams[22], affineParams[23]);
    *col = (*col + col3) / 2;
  }
}

extern void system1Sym(struct point * p, float * col, struct systemInfo * info) {
  float k = randomParam(0, 1);
  if (info->numFunctions == 2) {
    pickFrom2(p, col, info->functions, info->weights, info->affineParams, k, 0);
  } else if (info->numFunctions == 3) {
    pickFrom3(p, col, info->functions, info->weights, info->affineParams, k, 0);
  } else {
    pickFrom4(p, col, info->functions, info->weights, info->affineParams, k, 0);
  }
}

extern void system2Sym(struct point * p, float * col, struct systemInfo * info) {
  float k = randomParam(0, 1);
  if (k < 0.5) {
    // 180 degree rotation, no change in color to maintain symmetry
    *p = rotate(*p, (180 * M_PI) / 180);
  } else {
    if (info->numFunctions == 2) {
      pickFrom2(p, col, info->functions, info->weights, info->affineParams,
          k, 0.5);
    } else if (info->numFunctions == 3) {
      pickFrom3(p, col, info->functions, info->weights, info->affineParams,
          k, 0.5);
    } else {
      pickFrom4(p, col, info->functions, info->weights, info->affineParams,
          k, 0.5);
    }
  }
}

extern void system3Sym(struct point * p, float * col, struct systemInfo * info) {
  float k = randomParam(0, 1);
  if (k < 0.333) {
    // 120 degree rotation, no change in color to maintain symmetry
    *p = rotate(*p, (120 * M_PI) / 180);
  } else if (k < 0.667) {
    // 240 degree rotation, no change in color to maintain symmetry
    *p = rotate(*p, (240 * M_PI) / 180);
  } else {
    if (info->numFunctions == 2) {
      pickFrom2(p, col, info->functions, info->weights, info->affineParams,
          k, 0.667);
    } else if (info->numFunctions == 3) {
      pickFrom3(p, col, info->functions, info->weights, info->affineParams,
          k, 0.667);
    } else {
      pickFrom4(p, col, info->functions, info->weights, info->affineParams,
          k, 0.667);
    }
  }
}

