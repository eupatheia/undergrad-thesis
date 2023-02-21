#ifndef TRANSFORMS_H_
#define TRANSFORMS_H_

// keeps intermediate counts as we iterate the function
struct pix_counts {
  float sum;  // sum of accumulated color
  int alpha;  // defines final scaling
};

struct point {
  float x;
  float y;
};

extern struct point affine(struct point p, float a, float b, float c, float d,
    float e, float f);
extern struct point spherical(struct point p);

#endif
