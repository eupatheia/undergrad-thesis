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

// spherical variation 2
extern struct point spherical(struct point p) {
  float r = sqrt(pow(p.x, 2) + pow(p.y, 2));
  struct point next;
  next.x = p.x / pow(r, 2);
  next.y = p.y / pow(r, 2);
  return next;
}