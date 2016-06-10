/* PROGRAM gridlink1D

   --- gridlink np rmin rmax rcell z &ngrid &gridinit &gridlist
   --- Creates a 1D grid and places particles into it via a linked
   --- list.  Similar to gridlink.c, but in 1D.

   ---inputs---
      * np = number of particles
      * rmin,rmax = particles are located in a box running from 
                    (rmin,rmin,rmin) to (rmax,rmax,rmax).
      * rcell = size of a single grid cell 
      * z = array of particle coordinate that determines grid
   ---outputs---
      * ngrid = dimension of grid - computed from rmin,rmax,rcell
      * grid = 1D grid array where each cell contains the index of the 
                  first particle in that cell.
      * gridlist = array of length np containing linked list.
   -------------------------------------------------------------
      If cell (iz) contains N particles with indices j1,j2,j3,...,jN, 
      then: j1 = grid[iz], j2 = gridlist[j1], j3 = gridlist[j2],...,
      jN = gridlist[j<N-1>], and gridlist[jN] = -1.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>

#include "utils.h"
#include "cellarray.h"


#define NGRIDMAX 100000      /* maximum grid dimension */
#define MEMORY_INCREASE_FAC   1.5

void gridlink1D_with_struct(int np,double dmin,double dmax,double rcell,double *x1,double *y1,double *z1,double *dec, int *ngrid,cellarray **lattice)

{
  int nmesh,iz ;
  const double ddec = dmax-dmin;
  double inv_ddec = 1.0/ddec;
  cellarray *tmp;
  int expected_n,index;
  
  assert(ddec > 0.0);
  assert(rcell > 0.0);
  assert(MEMORY_INCREASE_FAC >= 1.0); 
  
  nmesh = (int)(ddec*BIN_REFINE_FACTOR/rcell) ;
  if(nmesh>NGRIDMAX) nmesh=NGRIDMAX ;
  *ngrid=nmesh ;

  expected_n=(int)((np/(double) nmesh)  *MEMORY_INCREASE_FAC);
  fprintf(stderr,"gridlink1D> Allocating %0.2g (MB) memory for the lattice, expected_n = %d nmesh = %d np=%d \n",(4*8+4)*expected_n*nmesh/(1024.*1024.),expected_n,nmesh,np);

  /*---Allocate-and-initialize-grid-arrays----------*/
  *lattice=(cellarray *) my_malloc(sizeof(cellarray),nmesh);
  for(int i=0;i<nmesh;i++) {
    tmp = &((*lattice)[i]);
    /* gettimeofday(&t0,NULL); */
    tmp->x     = my_malloc(sizeof(double),expected_n);
    tmp->y     = my_malloc(sizeof(double),expected_n);
    tmp->z     = my_malloc(sizeof(double),expected_n);
    tmp->dec    = my_malloc(sizeof(double),expected_n);
    tmp->index = my_malloc(sizeof(int),expected_n);
    tmp->nelements=0;
    tmp->nallocated=expected_n;
    /* gettimeofday(&t1,NULL); */
    /* malloc_time += ADD_DIFF_TIME(t0,t1); */
  }

  
  /*---Loop-over-particles-and-build-grid-arrays----*/
  for(int i=0;i<np;i++) {
    if(dec[i] >=dmin && dec[i] < dmax) {
      iz = (int)(nmesh*(dec[i]-dmin)*inv_ddec) ;  
      /* gettimeofday(&t0,NULL); */
      tmp = &((*lattice)[iz]);
      /* gettimeofday(&t1,NULL); */
    /* address_time += ADD_DIFF_TIME(t0,t1); */
      
    /* gettimeofday(&t00,NULL); */
      if(tmp->nelements == tmp->nallocated) {
	expected_n = tmp->nallocated*MEMORY_INCREASE_FAC;
	//fprintf(stderr,"allocated = %d new allocation = %d iz = %d \n",tmp->nallocated,expected_n,iz);
	
	/* gettimeofday(&t0,NULL); */
	tmp->x = my_realloc_in_function((void **) &(tmp->x) ,sizeof(double),expected_n,"lattice.x");
	tmp->y = my_realloc_in_function((void **) &(tmp->y) ,sizeof(double),expected_n,"lattice.y");
	tmp->z = my_realloc_in_function((void **) &(tmp->z) ,sizeof(double),expected_n,"lattice.z");
	tmp->dec= my_realloc_in_function((void **) &(tmp->dec),sizeof(double),expected_n,"lattice.dec");


	
	tmp->index=my_realloc_in_function((void **) &(tmp->index),sizeof(int),expected_n,"lattice.index");
	/* gettimeofday(&t1,NULL); */
	/* malloc_time += ADD_DIFF_TIME(t0,t1); */
	
	tmp->nallocated = expected_n;
      }
 
      
	/* gettimeofday(&t1,NULL); */
      /* if_time += ADD_DIFF_TIME(t00,t1); */
      
      /* gettimeofday(&t0,NULL); */
      index=tmp->nelements;
      tmp->x[index] = x1[i];
      tmp->y[index] = y1[i];
      tmp->z[index] = z1[i];
      tmp->dec[index] = dec[i];

      tmp->index[index] = i;
      tmp->nelements++;
      /* gettimeofday(&t1,NULL); */
    /* assignment_time += ADD_DIFF_TIME(t0,t1); */
    }
  }
 	
  
/* gettimeofday(&t1,NULL); */
  /* fprintf(stderr,"gridlink1D> Total = %0.3lf malloc_time = %0.3lf sec assignment_time = %0.3lf address_time = %0.3lf if_time = %0.3lf \n",ADD_DIFF_TIME(tbegin,t1),malloc_time,assignment_time,address_time,if_time); */
}


/* #define NGRIDMAX 100000      /\* maximum grid dimension *\/ */

/* void gridlink1D(int np,double rmin,double rmax,double rcell,double *z,int *ngrid,int **gridinit,int **gridlist) */
/* { */
/*   int nmesh,i,j,k,iz ; */

/*   nmesh = (int)((rmax-rmin)/rcell) ; */
/*   if(nmesh>NGRIDMAX) nmesh=NGRIDMAX ; */
/*   *ngrid=nmesh ; */

/* /\*---Allocate-and-initialize-grid-arrays----------*\/ */

/*   *gridlist=(int *)calloc(np,sizeof(int)) ; */
/*   *gridinit=(int *)calloc(nmesh,sizeof(int)) ; */
/*   for(i=0;i<np;i++) */
/*     (*gridlist)[i]=-1 ; */
/*   for(i=0;i<nmesh;i++) */
/*     (*gridinit)[i]=-1 ; */

/* /\*---Loop-over-particles-and-build-grid-arrays----*\/ */

/*   for(i=0;i<np;i++) */
/*     { */
/*       iz = (int)(nmesh*(z[i]-rmin)/(rmax-rmin)) ; */
/*       if(iz>=nmesh)  */
/* 	{ */
/* 	  /\* */
/* 	  fprintf(stderr,"gridlink> warning, particle %d at position z= %f\n",i,z[i]) ; */
/* 	  *\/ */
/* 	  iz-- ; */
/* 	} */
/*      (*gridlist)[i]=(*gridinit)[iz] ; */
/*      (*gridinit)[iz]=i ; */
/*     } */
/* } */






