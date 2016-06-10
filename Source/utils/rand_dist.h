#ifndef __RAND_DIST_H
#define __RAND_DIST_H

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "rng_gsl.h"
#include <gsl/gsl_randist.h>


int Average(double Nexp, void * rng );
int Poisson(double Nexp, void * rng );
int Binomial(double Nexp, void * rng);
int NegBinomial(double Nexp, void * rng );
inline int binary_search_probabilities(const double *prob,const int Npart,const double rand);
inline int binary_search_probabilities_float(const float *prob,const int Npart,const double rand);


#endif
