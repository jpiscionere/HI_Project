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

    int nthreads=16,chunk=CHUNKSIZE;

	int i, j, k;

	FILE *fp1, *fp2;
	char *galaxy_file, *group_file;

/* Input Properties of Galaxies */	
	
	char *Galaxy_Morphology;
	double *Galaxy_RA,
			*Galaxy_Dec,
			*Galaxy_Velocity,
			*Galaxy_Bj_mag,
			*Galaxy_R_mag,
			*Galaxy_I_mag;
			
			
/* Calculated Properties of Galaxies */
			
	int		N_Gal=0,
			Galaxy_Size=1E5;
	double	*Galaxy_Distance,
			*Galaxy_X,
			*Galaxy_Y,
			*Galaxy_Z,
			*Galaxy_HI_Mass;

/* Input Properties of Groups */

	int		N_Group=0,
			Group_Size=1E6;

	double *Group_RA,
			*Group_Dec,
			*Group_Velocity,
			*Group_Given_Distance,
			*Group_2_Radius,
			*Group_Vir_Radius,
			*Group_Ngroup;

/* Calculated Properties of Groups */

	double *Group_Distance,
			*Group_X,
			*Group_Y,
			*Group_Z,
			*Group_HI_Mass,
			*Group_Ltot;			
	
/*Pair Finding Quantities */
	
	int **Galaxies_In_Group,
  		*Belong_to_Group,
  		*Beong_to_Galaxy;
  	
  	double Min_Dec_Separation,
  			distance_SQR;
  		
  void gridlink1D_with_struct(int np,double dmin,double dmax,double rcell,double *x1,double *y1,double *z1,double *dec,int *ngrid,cellarray **lattice);

  struct timeval t0,t1;	

/* Randoms Counters */

	int nitems,
		nread,
		trash_d,
		flag;

	double trash_lf;
	
	char buffer[MAXBUFSIZE];


/* Read in Args */

	galaxy_file=argv[1];
	group_file=argv[2];


/* Read in Galaxy File */

	Galaxy_RA       = my_calloc(sizeof(*Galaxy_RA),Galaxy_Size);
	Galaxy_Dec      = my_calloc(sizeof(*Galaxy_Dec),Galaxy_Size);
	Galaxy_Velocity = my_calloc(sizeof(*Galaxy_Velocity),Galaxy_Size);
	Galaxy_Bj_mag   = my_calloc(sizeof(*Galaxy_Bj_mag),Galaxy_Size);
	Galaxy_R_mag    = my_calloc(sizeof(*Galaxy_R_mag),Galaxy_Size);
	Galaxy_I_mag    = my_calloc(sizeof(*Galaxy_I_mag),Galaxy_Size);
	Galaxy_Morphology      = my_calloc(sizeof(*Galaxy_Morphology),Galaxy_Size);
	
	  /////////////////////////////* [ READ IN THE GALAXY FILES AND CONVERT REDSHIFTS TO MPC ] *////////////////////////////////////
  	

	gettimeofday(&t0,NULL);
	fp1 = my_fopen(galaxy_file,"r") ;
	i=0;
	nitems=7;

	while(fgets(buffer,MAXBUFSIZE,fp1)!=NULL) {
    	nread=sscanf(buffer,"%lf %lf %lf %lf %lf %lf %c",
    		&Galaxy_RA[i],&Galaxy_Dec[i],&Galaxy_Velocity[i],
    		&Galaxy_Bj_mag[i],&Galaxy_R_mag[i], &Galaxy_I_mag[i],
    		&Galaxy_Morphology[i]);

		if (nread == nitems) {
      		i++;

			if(i==Galaxy_Size) {
	
				fprintf(stderr,"Increasing memory allocation for the spectroscopic sample\n");
				Galaxy_Size *= MEMORY_INCREASE_FAC;
				
				Galaxy_RA       = my_realloc(Galaxy_RA,sizeof(*Galaxy_RA),Galaxy_Size,"Galaxy_RA");
				Galaxy_Dec      = my_realloc(Galaxy_Dec,sizeof(*Galaxy_Dec),Galaxy_Size,"Galaxy_Dec");
				Galaxy_Velocity = my_realloc(Galaxy_Velocity,sizeof(*Galaxy_Velocity),Galaxy_Size,"Galaxy_Velocity");
				Galaxy_Bj_mag   = my_realloc(Galaxy_Bj_mag,sizeof(*Galaxy_Bj_mag),Galaxy_Size,"Galaxy_Bj_mag");
				Galaxy_R_mag    = my_realloc(Galaxy_R_mag,sizeof(*Galaxy_R_mag),Galaxy_Size,"Galaxy_R_mag");
				Galaxy_I_mag    = my_realloc(Galaxy_I_mag,sizeof(*Galaxy_I_mag),Galaxy_Size,"Galaxy_I_mag");
				Galaxy_Morphology      = my_realloc(Galaxy_Morphology,sizeof(*Galaxy_Morphology),Galaxy_Size,"Galaxy_Morphology");
	

      }
    } else {
      fprintf(stderr,"WARNING: In spectroscopic sample line %d did not contain %d elements...skipping line\n",i,nitems);
    }
  }
  N_Gal=i;
  fclose(fp1);
  gettimeofday(&t1,NULL);


/* Calculate Galaxy Properties */



/* Read in Group File */


/* Calculate Group Properties */










  





return 0;
}