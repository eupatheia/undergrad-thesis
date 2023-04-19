#ifndef MAT_H_
#define MAT_H_

#include "vec.h"

// representation of a 4x4 matrix (as 4 column vectors)
//    ax bx cx dx
//    ay by cy dy
//    az bz cz dz
//    aw bw cw dw
struct mat {
  struct vec a;
  struct vec b;
  struct vec c;
  struct vec d;
};

void printMat(const struct mat m);

struct mat matMul(const struct mat m, const struct mat n);
struct vec trans(const struct mat m, const struct vec v);
struct mat transpose(const struct mat m);

struct mat getViewMat(const struct vec eyePos, const struct vec lookPos,
    const struct vec upDir);

#endif