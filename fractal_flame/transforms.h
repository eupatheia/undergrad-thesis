#ifndef TRANSFORMS_H_
#define TRANSFORMS_H_

// function pointer type
typedef struct point (*transform_ptr) (struct point, float, float,
    float, float, float, float);

// keeps intermediate counts as we iterate the function
struct pix_counts {
  float color;  // blended color after all iterations, in range [0,1]
  int alpha;  // defines final scaling
};

struct point {
  float x;
  float y;
};

struct systemInfo {
  int numFunctions;
  int symmetry;
  transform_ptr * functions;
  float * weights;
  float * affineParams;
};

// ============== BASIC TRANSFORMS & HELPERS ==============
extern struct point affine(struct point p, float a, float b, float c, float d,
    float e, float f);
extern struct point rotate(struct point p, float angle);
extern float randomParam(float start, float end);

// ============== VARIATIONS (with Draves numbers) ==============
extern struct point linear(struct point p, float a, float b, float c, float d,
    float e, float f);  // var 0
extern struct point sinusoidal(struct point p, float a, float b, float c,
    float d, float e, float f);  // var 1
extern struct point spherical(struct point p, float a, float b, float c,
    float d, float e, float f);  // var 2
extern struct point swirl(struct point p, float a, float b, float c, float d,
    float e, float f);  // var 3
extern struct point horseshoe(struct point p, float a, float b, float c, float d,
    float e, float f);  // var 4
extern struct point polar(struct point p, float a, float b, float c, float d,
    float e, float f);  // var 5
extern struct point handkerchief(struct point p, float a, float b, float c,
    float d, float e, float f);  // var 6
extern struct point heart(struct point p, float a, float b, float c,
    float d, float e, float f);  // var 7
extern struct point disk(struct point p, float a, float b, float c,
    float d, float e, float f);  // var 8
extern struct point spiral(struct point p, float a, float b, float c, float d,
    float e, float f);  // var 9
extern struct point hyperbolic(struct point p, float a, float b, float c,
    float d, float e, float f);  // var 10
extern struct point diamond(struct point p, float a, float b, float c,
    float d, float e, float f);  // var 11
extern struct point ex(struct point p, float a, float b, float c, float d,
    float e, float f);  // var 12
extern struct point julia(struct point p, float a, float b, float c, float d,
    float e, float f);  // var 13
extern struct point fisheye(struct point p, float a, float b, float c, float d,
    float e, float f);  // var 16
extern struct point exponential(struct point p, float a, float b, float c,
    float d, float e, float f);  // var 18
extern struct point eyefish(struct point p, float a, float b, float c, float d,
    float e, float f);  // var 27

// ============== SYSTEMS ==============
extern void systemSphericalSpherical(struct point * p, float * c);
extern void systemLinearSwirlSpiral(struct point * p, float * c);
extern void systemLinearSwirlSpiral_3SYM(struct point * p, float * c);

// helper to choose among functions based on random val k and 
// booster value (leftover weight reserved for rotation functions)
extern void pickFrom2(struct point * p, float * col, transform_ptr * functions,
    float * weights, float * affineParams, float k, float booster);
extern void pickFrom3(struct point * p, float * col, transform_ptr * functions,
    float * weights, float * affineParams, float k, float booster);
extern void pickFrom4(struct point * p, float * col, transform_ptr * functions,
    float * weights, float * affineParams, float k, float booster);

// random function systems (1, 2, or 3-way symmetry)
extern void system1Sym(struct point * p, float * col, struct systemInfo * info);
extern void system2Sym(struct point * p, float * col, struct systemInfo * info);
extern void system3Sym(struct point * p, float * col, struct systemInfo * info);



#endif
