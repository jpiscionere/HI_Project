#ifndef __RNG_GSL_H
#define __RNG_GSL_H
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "utils.h"

void *rng_create( unsigned long int seed );
void *rng_create_per_thread(void);
void rng_free( void * state );
double rng_uniform( void * state );
double rng_gaussian( void * state , double sigma );


#endif
