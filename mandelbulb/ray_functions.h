// helper functions (SDFs, intersectors, raymarchers, normals, etc.)
// for mandelbulb.c

#ifndef RAY_FUNCTIONS_H_
#define RAY_FUNCTIONS_H_

float clampFloat(float num, float start, float end);
int clampInt(int num, int start, int end);
float density(struct vec localPos);
struct vec puffcolor(struct vec localPos);
int findFace(struct vec pos, float left, float right, float bottom, float top,
    float near, float far);
struct vec f(struct vec origPos, struct vec v);
struct vec optf(struct vec origPos, struct vec v);

float sdScene(struct vec pos, int maxIterations);
float sdBox(struct vec pos, struct vec ray, float xmin, float xmax, float ymin,
    float ymax, float zmin, float zmax);
struct vec sphereIntersect(struct vec pos, struct vec ray, struct vec center,
    float radius);

struct vec calcNormal(struct vec pos, float epsilon, int maxIterations);

int DAraymarch(struct vec pos, struct vec ray, struct vec * hitPos,
    float * shadowColor, int maxSteps, float hitRange, int maxIterations,
    float softness);
int STEPraymarch(struct vec pos, struct vec ray, struct vec * hitPos,
    float maxMarchDistance, float stepSize, int maxIterations);

struct vec farBulbIntersect(struct vec origPos, struct vec ray,
    struct vec hitPos, int maxIterations, int maxSteps, float hitRange,
    float softness);

#endif