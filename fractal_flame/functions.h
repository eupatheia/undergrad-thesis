#ifndef FUNCTIONS_H_
#define FUNCTIONS_H_

// keeps intermediate counts as we iterate the function
struct pix_counts {
  // summed over all iterations, will be scaled to get RBG color at end
  float countR;
  float countG;
  float countB;
  int countA;  // alpha channel, defines final scaling
};

struct point {
  float x;
  float y;
};

extern struct point affine(struct point p, float a, float b, float c, float d,
    float e, float f);
extern struct point spherical(struct point p);

#endif
