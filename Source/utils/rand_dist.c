/* These are modified routines from the original to be a little more efficient
 * as well as use a "custom" wrapper for the uniform distribution from GSL */

/* requires "rng_gsl.c" to also be included */

// XXX these need to be tested!

/*-Average-distribution------------------------------*/

#include "rand_dist.h"


int binary_search_probabilities(const double *prob,const int Npart,const double rand)
{
  const int threshold_for_linear_search=6;
  int left, right, middle;
  left = 0;
  right = Npart-1;
  
  /* while(left <= right) { */
  while((right-left) > threshold_for_linear_search) {
    middle = (left + right) >> 1;//divide by 2                                                                                                                                                                                                 
    if(rand <= prob[middle]) {
      right = middle-1;
    } else {
      left = middle + 1;
    }

    /* if((right-left) < threshold_for_linear_search) { */
    /*   break; */
    /* } */
  
  }
  
  /* fprintf(stderr,"left = %d right=%d \n",left,right); */
  for(middle=left;middle<=right;middle++)
    if(rand <= prob[middle])
      break;
  
  return middle;
}


int binary_search_probabilities_float(const float *prob,const int Npart,const double rand)
{
  const int threshold_for_linear_search=6;
  int left, right, middle;
  left = 0;
  right = Npart-1;
  
  /* while(left <= right) { */
  while((right-left) > threshold_for_linear_search) {
    middle = (left + right) >> 1;//divide by 2                                                                                                                                                                                                 
    if(rand <= prob[middle]) {
      right = middle-1;
    } else {
      left = middle + 1;
    }

    /* if((right-left) < threshold_for_linear_search) { */
    /*   break; */
    /* } */
  
  }
  
  /* fprintf(stderr,"left = %d right=%d \n",left,right); */
  for(middle=left;middle<=right;middle++)
    if(rand <= prob[middle])
      break;
  
  return middle;
}




int Average(double Nexp, void * rng )
{
  int Nact ;
  double rand ;

  rand = rng_uniform(rng);

  /* if(rand<=(Nexp-(int)(Nexp))) Nact = (int) (Nexp +1) ; */
  /* else Nact = (int)(Nexp) ; */

  if(rand<=(Nexp-(int)(Nexp))) Nact = (int) ceil(Nexp);
  else Nact = (int) floor(Nexp) ;

  
  return Nact ;
}

/*-Poisson-distribution------------------------------*/
int Poisson(double Nexp, void * rng )
{

  return gsl_ran_poisson( rng, Nexp);
  
  /* int i,Nmax,Nact=0 ; */
  /* double sigma,P,rand ; */
  /* double sum; */

  /* sigma = sqrt(Nexp) ; */
  /* if(Nexp>=0.6) Nmax = (int)(10*sigma+Nexp) ; */
  /* else Nmax = 8; */

  /* P = exp(-Nexp) ; */
  /* sum = P ; */

  /* rand = rng_uniform(rng); */

  /* i = 0; */
  /* while(i < Nmax) */
  /*   { */
  /*     if(rand<=sum) */
  /*       { */
  /*         Nact = i ; */
  /*         break ; */
  /*       } */
  /*     i++ ; */
  /*     P *= Nexp/(double)i; */
  /*     sum += P; */
  /*   } */

  /* return Nact ; */
}

/*-Binomial-distribution-----------------------------*/
int Binomial(double Nexp, void * rng)
{
  int i,Nmax,Nact=0,r ;
  double x,p,q,Prob;//,var ;
  double rand ;
  double sum ;

  r=(int)(2*Nexp+1) ;
  p=Nexp/(double)r ;
  q=1-p ;
  //var=sqrt(r*p*q) ;
  if(Nexp>=4) Nmax = r;
  else Nmax=8 ;

  Prob = pow(q,(double)r) ;
  sum = Prob ;

  rand = rng_uniform(rng);

  i = 0;
  while(i < Nmax)
    {
      if(rand<=sum)
        {
          Nact = i ;
          break ;
        }
      i++ ;
      x = (double)i ;
      Prob *= ((r-x+1)/x)*p/q ;
      sum += Prob;
    }

  return Nact ;
}

/*-Negative-Binomial-distribution--------------------*/
int NegBinomial(double Nexp, void * rng )
{
  int i,Nmax,Nact=0,r ;
  double x,p,q,Prob,var ;
  double rand ;
  double sum ;

  r=(int)Nexp+1 ;
  p=1/(1+Nexp/(double)r) ;
  q=1-p ;
  var=sqrt(r*q/(p*p)) ;

  if(Nexp>=3) Nmax = (int)(10*var+Nexp) ;
  else Nmax = 30;

  Prob = pow(p,(double)r) ;
  sum = Prob ;

  rand = rng_uniform( rng );

  i = 0;
  while(i < Nmax)
    {
      if(rand<=sum)
        {
          Nact = i ;
          break ;
        }
      i++ ;
      x = (double)i ;
      Prob *= ((r+x-1)/x)*pow(q,x)/pow(q,x-1) ;
      if(i>r)
        {
          if(Prob>1e-20) Prob=Prob ;
          else Prob=0 ;
        }
      sum += Prob ;
    }

  return Nact ;
}


