/* this just simplifies using the GSL RNG */
#include "rng_gsl.h"

/* initialize GSL_RNG and return a "state" variable */
void *rng_create( unsigned long int seed )
{
    gsl_rng * rng;

    /* this sets the type, MT19937 is a good RNG */
    rng = gsl_rng_alloc( gsl_rng_mt19937 );
#ifdef DEBUG
    fprintf(stderr, "initializing random number generator using seed=%lu\n", seed);
#endif    
    gsl_rng_set(rng, seed);

    return (void *)rng;
}

void *rng_create_per_thread(void)
{
    gsl_rng * rng;
    unsigned long seed;
    const char fname[]="/dev/urandom";
    FILE *fp=my_fopen(fname, "r");
    my_fread(&seed, sizeof(seed),1,fp);
    fclose(fp);

    /* this sets the type, MT19937 is a good RNG */
    rng = gsl_rng_alloc( gsl_rng_mt19937 );
#ifdef DEBUG    
    fprintf(stderr, "initializing random number generator using seed=%lu\n", seed);
#endif    
    gsl_rng_set(rng, seed);

    return (void *)rng;
}


/* clean up to deallocate RNG structure */
void rng_free( void * state )
{
    gsl_rng_free( (gsl_rng *) state );
}

/* a wrapper to yield a uniform deviate */
double rng_uniform( void * state )
{
    return gsl_rng_uniform( (gsl_rng *) state );
}

/* a wrapper to yield a gaussian deviate */
double rng_gaussian( void * state , double sigma )
{
    return gsl_ran_gaussian( (gsl_rng *) state , sigma );
}
