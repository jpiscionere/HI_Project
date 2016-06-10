#ifndef __REDSHIFT_DISTORTIONS_H_
#define __REDSHIFT_DISTORTIONS_H_
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_interp.h>

#include "utils.h"
#include "cosmology_params.h"
#include "set_cosmo_dist.h"


#define sqr(x) ((x)*(x))
#define cube(x) ((x)*(x)*(x))


#define AXIS_OFFSET 60.0 
/* #define ZMAX 2.0 */
/* #define Ho   100.0  */

void pp_zdistortions(float *z, float *vz,float Lbox, const int np, const int zspace,const float zmedian, const int lasdamas_cosmology);

#endif
