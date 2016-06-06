#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <gsl/gsl_integration.h>
#include <string.h>

#include "cellarray.h"

#include "utils.h"
#include "progressbar.h"
#include <omp.h>
#define CHUNKSIZE   10

#define MAXLEN          (5000)
#define MAXBUFSIZE      (10000)

//Change these values to the correct double-precision ones
#define PI              (3.141592)
#define DEG_TO_RAD      (0.01745328888)
#define RAD_TO_DEG      (57.2957914331)


#define SQR(x) ((x)*(x))
#define SQRT(x) (pow(x,0.5))
#define OMEGA_M (0.3)
#define OMEGA_L (0.7)
#define W_INDEX (-1.0)
#define H_0 (100)
#define SPEED_OF_LIGHT (299792)
#define LITTLE_H (1)


#define MEMORY_INCREASE_FAC                               (1.2)


int main(int argc, char *argv[])
{

int i, j, k;

  FILE *fp1, *fp2;
  char *galaxy_file_1, *galaxy_file_2;








return 0;
}