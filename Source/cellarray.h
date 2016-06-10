#ifndef __CELLARRAY_H
#define __CELLARRAY_H

#include <stdio.h>
#include <stdlib.h>

#define BIN_REFINE_FACTOR     4

#define REGISTER_WIDTH 256  //cpu supports avx instructions
#define NVECF  8  //8 floats per ymm register
#define NVECD  4  //4 doubles per ymm register

#define ALIGNMENT 32

//cellarray does not "need" to be aligned but the individual elements do.
typedef struct{
  double *x  __attribute__((aligned(ALIGNMENT)));
  double *y  __attribute__((aligned(ALIGNMENT)));
  double *z  __attribute__((aligned(ALIGNMENT)));
  double *dec __attribute__((aligned(ALIGNMENT)));
  double *mag __attribute__((aligned(ALIGNMENT)));
  int *index __attribute__((aligned(ALIGNMENT)));
  int nelements;
  int nallocated;
}cellarray __attribute__((aligned(ALIGNMENT)));

#endif
