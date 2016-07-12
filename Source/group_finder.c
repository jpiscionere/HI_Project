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
#define H_0 (100.)
#define SPEED_OF_LIGHT (299792.)
#define LITTLE_H (1.)

#define 	GSL_MAX(a, b)   ((a) > (b) ? (a) : (b))
#define 	GSL_MIN(a, b)   ((a) < (b) ? (a) : (b))


#define ADD_DIFF_TIME(tstart,tend)                        ((tend.tv_sec - tstart.tv_sec) + 1e-6*(tend.tv_usec-tstart.tv_usec)) 
#define MEMORY_INCREASE_FAC                               (1.2)


int main(int argc, char *argv[])
{

    int nthreads=4,chunk=CHUNKSIZE;

	int i, j, k;

	FILE *fp1, *fp2;
	char *galaxy_file, *group_file;

/* Input Properties of Galaxies */	
	
	char *Galaxy_Morphology;
	int *Galaxy_ID,
		*Belong_To_Group;
	double *Galaxy_RA,
			*Galaxy_Dec,
			*Galaxy_Velocity,
			*Galaxy_Bj_mag,
			*Galaxy_R_mag,
			*Galaxy_I_mag,
			*Galaxy_Sint;
			
			
/* Calculated Properties of Galaxies */
			
	int		N_Galaxy=0,
			Galaxy_Size=1E4;
	double	*Galaxy_Distance,
			*Galaxy_X,
			*Galaxy_Y,
			*Galaxy_Z,
			*Galaxy_HI_Mass;

/* Input Properties of Groups */

	int		N_Group=0,
			Group_Size=1E4,
			*Group_Ngroup;

	double *Group_RA,
			*Group_Dec,
			*Group_Velocity,
			*Group_Given_Distance,
			*Group_2_Radius,
			*Group_Radius_SQR;

/* Calculated Properties of Groups */

	double *Group_Distance,
			*Group_X,
			*Group_Y,
			*Group_Z,
			*Group_HI_Mass,
			*Group_Ltot;			
	
/*Pair Finding Quantities */
	
	int **Galaxies_In_Group,
		Max_Group_Members=20,
  		*Belong_To_Galaxy;
  	
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


/* Allocation For Galaxies */

	Galaxy_ID       = my_calloc(sizeof(*Galaxy_ID),Galaxy_Size);
	Galaxy_RA       = my_calloc(sizeof(*Galaxy_RA),Galaxy_Size);
	Galaxy_Dec      = my_calloc(sizeof(*Galaxy_Dec),Galaxy_Size);
	Galaxy_Velocity = my_calloc(sizeof(*Galaxy_Velocity),Galaxy_Size);
	Galaxy_Bj_mag   = my_calloc(sizeof(*Galaxy_Bj_mag),Galaxy_Size);
	Galaxy_R_mag    = my_calloc(sizeof(*Galaxy_R_mag),Galaxy_Size);
	Galaxy_I_mag    = my_calloc(sizeof(*Galaxy_I_mag),Galaxy_Size);
	Galaxy_Sint		= my_calloc(sizeof(*Galaxy_Sint),Galaxy_Size);
	Galaxy_Morphology      = my_calloc(sizeof(*Galaxy_Morphology),Galaxy_Size);


	Galaxy_X      = my_calloc(sizeof(*Galaxy_X),Galaxy_Size);
	Galaxy_Y      = my_calloc(sizeof(*Galaxy_Y),Galaxy_Size);
	Galaxy_Z      = my_calloc(sizeof(*Galaxy_Z),Galaxy_Size);
	Galaxy_Distance      = my_calloc(sizeof(*Galaxy_Distance),Galaxy_Size);
	Galaxy_HI_Mass      = my_calloc(sizeof(*Galaxy_HI_Mass),Galaxy_Size);

/* Allocation for Groups */

	Group_RA	= my_calloc(sizeof(*Group_RA),Group_Size);
	Group_Dec	= my_calloc(sizeof(*Group_Dec),Group_Size);
	Group_Velocity	= my_calloc(sizeof(*Group_Velocity),Group_Size);
	Group_Given_Distance	= my_calloc(sizeof(*Group_Given_Distance),Group_Size);
	Group_Ngroup	= my_calloc(sizeof(*Group_Ngroup),Group_Size);
	Group_2_Radius	= my_calloc(sizeof(*Group_2_Radius),Group_Size);

	Group_Radius_SQR	= my_calloc(sizeof(*Group_Radius_SQR),Group_Size);	
	Group_Distance	= my_calloc(sizeof(*Group_Distance),Group_Size);
	Group_X	= my_calloc(sizeof(*Group_X),Group_Size);
	Group_Y	= my_calloc(sizeof(*Group_Y),Group_Size);
	Group_Z	= my_calloc(sizeof(*Group_Z),Group_Size);
	Group_HI_Mass	= my_calloc(sizeof(*Group_HI_Mass),Group_Size);
	Group_Ltot	= my_calloc(sizeof(*Group_Ltot),Group_Size);
	
	
	  /////////////////////////////* [ READ IN THE GALAXY FILE ] *////////////////////////////////////
  	

	gettimeofday(&t0,NULL);
	fp1 = my_fopen(galaxy_file,"r") ;
	i=0;
	nitems=9;

	while(fgets(buffer,MAXBUFSIZE,fp1)!=NULL) {
    	nread=sscanf(buffer,"%d %lf %lf %lf %lf %lf %lf %lf %c",
    		&Galaxy_ID[i],&Galaxy_RA[i],&Galaxy_Dec[i],&Galaxy_Velocity[i],
    		&Galaxy_Bj_mag[i],&Galaxy_R_mag[i], &Galaxy_I_mag[i],
    		&Galaxy_Sint[i],&Galaxy_Morphology[i]);


		if (nread == nitems) {
      		Galaxy_Distance[i]=Galaxy_Velocity[i] / H_0;
      		Galaxy_X[i]= Galaxy_Distance[i] * sin((90-Galaxy_Dec[i])*DEG_TO_RAD)*cos(Galaxy_RA[i]*DEG_TO_RAD);
      		Galaxy_Y[i]= Galaxy_Distance[i] * sin((90-Galaxy_Dec[i])*DEG_TO_RAD)*sin(Galaxy_RA[i]*DEG_TO_RAD);
      		Galaxy_Z[i]= Galaxy_Distance[i] * cos((90-Galaxy_Dec[i])*DEG_TO_RAD);
      		Galaxy_HI_Mass[i] = 2.365 * 1E5 * Galaxy_Distance[i] * Galaxy_Distance[i] * Galaxy_Sint[i];
      		i++;

			if(i==Galaxy_Size) {
	
				fprintf(stderr,"Increasing memory allocation for the galaxy sample\n");
				Galaxy_Size *= MEMORY_INCREASE_FAC;
				
				Galaxy_RA       = my_realloc(Galaxy_RA,sizeof(*Galaxy_RA),Galaxy_Size,"Galaxy_RA");
				Galaxy_Dec      = my_realloc(Galaxy_Dec,sizeof(*Galaxy_Dec),Galaxy_Size,"Galaxy_Dec");
				Galaxy_Velocity = my_realloc(Galaxy_Velocity,sizeof(*Galaxy_Velocity),Galaxy_Size,"Galaxy_Velocity");
				Galaxy_Bj_mag   = my_realloc(Galaxy_Bj_mag,sizeof(*Galaxy_Bj_mag),Galaxy_Size,"Galaxy_Bj_mag");
				Galaxy_R_mag    = my_realloc(Galaxy_R_mag,sizeof(*Galaxy_R_mag),Galaxy_Size,"Galaxy_R_mag");
				Galaxy_I_mag    = my_realloc(Galaxy_I_mag,sizeof(*Galaxy_I_mag),Galaxy_Size,"Galaxy_I_mag");
				Galaxy_Morphology      = my_realloc(Galaxy_Morphology,sizeof(*Galaxy_Morphology),Galaxy_Size,"Galaxy_Morphology");

				Galaxy_Distance       = my_realloc(Galaxy_Distance,sizeof(*Galaxy_Distance),Galaxy_Size,"Galaxy_Distance");
				Galaxy_X       = my_realloc(Galaxy_X,sizeof(*Galaxy_X),Galaxy_Size,"Galaxy_X");
				Galaxy_Y       = my_realloc(Galaxy_Y,sizeof(*Galaxy_Y),Galaxy_Size,"Galaxy_Y");
				Galaxy_Z       = my_realloc(Galaxy_Z,sizeof(*Galaxy_Z),Galaxy_Size,"Galaxy_Z");
				Galaxy_HI_Mass       = my_realloc(Galaxy_HI_Mass,sizeof(*Galaxy_HI_Mass),Galaxy_Size,"Galaxy_HI_Mass");


      }
    } else {
      fprintf(stderr,"WARNING: In spectroscopic sample line %d did not contain %d elements...skipping line\n",i,nitems);
    }
  }
  N_Galaxy=i;
  fclose(fp1);
  gettimeofday(&t1,NULL);


  fprintf(stderr,"Group Finder > There are %d Galaxies. Time taken = %6.2lf sec\n",N_Galaxy,ADD_DIFF_TIME(t0,t1));	


	  /////////////////////////////* [ READ IN THE GROUP FILE ] *////////////////////////////////////


	gettimeofday(&t0,NULL);
	fp2 = my_fopen(group_file,"r") ;
	i=0;
	nitems=6;
	double	Max_Radius=0,Distance_to_Nearest_Group=20000;



	while(fgets(buffer,MAXBUFSIZE,fp2)!=NULL) {
		nread=sscanf(buffer,"%lf %lf %lf %lf %lf %d",
			&Group_RA[i],&Group_Dec[i],&Group_Velocity[i],&Group_Given_Distance[i],&Group_2_Radius[i],&Group_Ngroup[i]);
	
		if (nread == nitems) {
		

      		Group_Distance[i]=Group_Velocity[i] / H_0;
      		Group_X[i]= Group_Distance[i] * sin((90-Group_Dec[i])*DEG_TO_RAD)*cos(Group_RA[i]*DEG_TO_RAD);
      		Group_Y[i]= Group_Distance[i] * sin((90-Group_Dec[i])*DEG_TO_RAD)*sin(Group_RA[i]*DEG_TO_RAD);
      		Group_Z[i]= Group_Distance[i] * cos((90-Group_Dec[i])*DEG_TO_RAD);
      		Group_Radius_SQR[i]=Group_2_Radius[i]*Group_2_Radius[i];
      		
      		
      		if(Group_Distance[i] < Distance_to_Nearest_Group)
      			Distance_to_Nearest_Group = Group_Distance[i];

			if(Group_2_Radius[i] > Max_Radius)
				Max_Radius=Group_2_Radius[i];      		
      		
      		i++;

			if(i==Group_Size) {
	
				fprintf(stderr,"Increasing memory allocation for the group sample\n");
				Group_Size *= MEMORY_INCREASE_FAC;
				
				Group_RA       = my_realloc(Group_RA,sizeof(*Group_RA),Group_Size,"Group_RA");
				Group_Dec      = my_realloc(Group_Dec,sizeof(*Group_Dec),Group_Size,"Group_Dec");
				Group_Velocity = my_realloc(Group_Velocity,sizeof(*Group_Velocity),Group_Size,"Group_Velocity");
				Group_Given_Distance = my_realloc(Group_Given_Distance,sizeof(*Group_Given_Distance),Group_Size,"Group_Given_Distance");
				Group_2_Radius = my_realloc(Group_2_Radius,sizeof(*Group_2_Radius),Group_Size,"Group_2_Radius");
				Group_Radius_SQR = my_realloc(Group_Radius_SQR,sizeof(*Group_Radius_SQR),Group_Size,"Group_Radius_SQR");
				Group_Ngroup = my_realloc(Group_Ngroup,sizeof(*Group_Ngroup),Group_Size,"Group_Ngroup");

				Group_Distance       = my_realloc(Group_Distance,sizeof(*Group_Distance),Group_Size,"Group_Distance");
				Group_X       = my_realloc(Group_X,sizeof(*Group_X),Group_Size,"Group_X");
				Group_Y       = my_realloc(Group_Y,sizeof(*Group_Y),Group_Size,"Group_Y");
				Group_Z       = my_realloc(Group_Z,sizeof(*Group_Z),Group_Size,"Group_Z");

      		}
    	} else {
      		fprintf(stderr,"WARNING: In spectroscopic sample line %d did not contain %d elements...skipping line\n",i,nitems);
    	}
	}

	N_Group=i;
	fclose(fp2);
	gettimeofday(&t1,NULL);
  	fprintf(stderr,"Group Finder > There are %d Groups. Time taken = %6.2lf sec\n",N_Group,ADD_DIFF_TIME(t0,t1));	
	fprintf(stderr,"The Maximum Radius = %lf\n",Max_Radius);
	

// Grid Link Galaxies 

	int interrupted=0;
	int ngrid;/* *gridinit1D,*gridlist1D ; */
	double dmin=-90,dmax=0;//min/max dec
	double inv_dmax_diff = 1.0/(dmax-dmin);
	cellarray *lattice;
	cellarray *cellstruct __attribute__((aligned(ALIGNMENT)));
	init_my_progressbar(N_Group,&interrupted);
    omp_set_num_threads(nthreads);

  
  	int xx=0;
  	

  	double Radius_Cut=2.* Max_Radius;
  	fprintf(stderr,"Radius Cut = %lf\n",Radius_Cut);
  	Galaxies_In_Group    = my_malloc(sizeof(int *),N_Group);
  	
  	
  	for(i=0;i<N_Group;i++)
    	Galaxies_In_Group[i] = my_calloc(sizeof(int),Max_Group_Members);
	
	Belong_To_Group = my_calloc(sizeof(*Belong_To_Group),N_Galaxy);
	Group_HI_Mass	= my_calloc(sizeof(*Group_HI_Mass),N_Group);
	int i_flag=0,group_counter=0;
/*
#pragma omp parallel default(none) shared(interrupted,N_Galaxy,N_Group,Group_Distance,Galaxy_Distance,Radius_Cut,Group_X,Group_Y,Group_Z,Galaxy_X,Galaxy_Y,Galaxy_Z,Galaxies_In_Group,Group_HI_Mass,Belong_To_Group)
{
	int tid = omp_get_thread_num();
	#pragma omp for schedule(dynamic,chunk)
*/
	for(i=0;i<N_Group;i++){
		init_my_progressbar(N_Group,&interrupted);
		if(Group_Velocity[i]<=1.0E5 && Group_Ngroup[i] > 4){
			for(j=0;j<N_Galaxy;j++){
				if(fabs(Group_Distance[i] - Galaxy_Distance[j]) < Radius_Cut){
					distance_SQR=(Group_X[i] - Galaxy_X[j])*(Group_X[i] - Galaxy_X[j]) +
								(Group_Y[i] - Galaxy_Y[j])*(Group_Y[i] - Galaxy_Y[j]) +
								(Group_Z[i] - Galaxy_Z[j])*(Group_Z[i] - Galaxy_Z[j]);
							
					if(distance_SQR < Group_Radius_SQR[i]){
						Galaxies_In_Group[i][0]++;
						Galaxies_In_Group[i][Galaxies_In_Group[i][0]] = j;
						Group_HI_Mass[i]+=Galaxy_HI_Mass[j];
						Belong_To_Group[j]++;
						i_flag=1;
					}
				
				} 
			}
		}		
	}
//}	

	for(i=0;i<N_Group;i++)
		if(Galaxies_In_Group[i][0] > 0)
			group_counter++;
	
	int galaxy_counter=0;
	int belong_to_multiple_groups=0;
	int max_group=0,max_galaxy=0;
	
	for(i=0;i < N_Galaxy;i++){
		galaxy_counter += Belong_To_Group[i];
	
		if(Belong_To_Group[i] > 1)
			belong_to_multiple_groups += Belong_To_Group[i];
		if(Belong_To_Group[i] > max_group)
			max_group=Belong_To_Group[i];		
	}
	
	fprintf(stderr,"There are %d Groups with Galaxies Associated with them\n",group_counter);
	fprintf(stderr,"%d Galaxies are in Groups (%lf of total). %d Galaxies are in multiple groups\n",
					galaxy_counter,galaxy_counter*1.0/N_Galaxy,belong_to_multiple_groups);
	fprintf(stderr,"The Maximum number of Groups a Galaxy Belongs to is %d\n",max_group);
	 
     finish_myprogressbar(&interrupted);


//	for(i=0;i<N_Galaxy;i++){
//		fprintf(stdout,"%lf %lf %lf %lf %lf %d",Galaxy_ID[i],Galaxy_RA[i],Galaxy_Dec[i],Galaxy_Velocity[i],
//   		Galaxy_Bj_mag[i],Galaxy_R_mag[i], Galaxy_I_mag[i],
//    		Galaxy_Sint[i],Galaxy_Morphology[i]);
	
//	}



return 0;
}