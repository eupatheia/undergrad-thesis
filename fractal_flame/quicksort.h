#ifndef QUICKSORT_H_
#define QUICKSORT_H_

#include <stdio.h>
#include <stdlib.h>
#include "transforms.h"

extern void quicksort(struct specimen * arr, int left, int right);
extern void hybridSort(struct specimen * arr, int left, int right);

#endif