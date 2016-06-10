/* PROGRAM GETMULT3

   --- getmult3 volume [Nmax] < groupfile > multfile
   --- Computes group multiplicity function
       Verison 3: reads output of fof4.c (.pnm file)

      * volume : volume of group sample (in h^3/Mpc^3)
      * [Nmax] : maximum N for which ngrp is evaluated (default=200) 
      * < groupfile : ascii file containing <ID N RA DEC cz czdisp Mgtot Mrtot comp redge> (fof4.c)
      * > output : ascii file containing <N Ngrps log(ngrp1) log(ngrp2) indx>
                   Ngrps : number of groups with =N
                   ngrp1 : space density of groups with >=N
                   ngrp2 : space density of groups with > N
		   indx  : (0,1,2) if (real, interpolated, extrapolated)
## Compile string: icc -std=c99 -xhost utils.c getmult3.c -o getmult3 -Wall -Wextra -Wshadow -g -I .
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "sglib.h"
#include "utils.h"

int main(int argc,char **argv)
{
  int Ngroups,*Ng=NULL,nmin;
  double volume;
  int i,j ;
  /*---Arguments-------------------------*/
  int Nmax;


  //For file i/o
  const int MAXBUFSIZE=10000;
  char *fname=NULL;
  char buffer[MAXBUFSIZE];
  const char comment='#';
  const int nitems = 1;
  int nread;
  FILE *fp = NULL;
  
  /*-------------------------------------*/
  int N,*indx,Nlast,*Ngrps ;
  double *logngrp1,*logngrp2 ;
  double X1,X2,Y1,Y2,Z1,Z2 ;
  double inv_volume,half_log10_inv_volume;
  int first_index_with_nmin=0;
  if(argc != 4) {
    fprintf(stderr,"ERROR: Usage - `%s' <total volume> <GalaxyFile> <nmin>",argv[0]);
    exit(EXIT_FAILURE);
  }
  
  volume=atof(argv[1]);
  fname=argv[2];
  nmin =atoi(argv[3]);
  
  Ngroups = getnumlines(fname,comment);
  Ng = my_malloc(sizeof(*Ng),Ngroups);
  fp = my_fopen(fname,"r");
  i=0;
  while(1) {
    if(fgets(buffer, MAXBUFSIZE,fp)!=NULL) {
      if(buffer[0] !=comment) {
	nread = sscanf(buffer,"%*s %*d %d ",&Ng[i]);
	assert(nread == nitems);
	i++;
      }
    } else
      break;
  }
  assert(i==Ngroups);
  fclose(fp);


  /*---Read-Arguments---------------------------------------------------*/
  assert(volume > 0.0);
  inv_volume=1.0/volume;
  half_log10_inv_volume = log10(0.5*inv_volume);  

  /*---Sort-groups-by-N-------------------------------------------------*/
  SGLIB_ARRAY_SINGLE_HEAP_SORT(int, Ng, Ngroups, SGLIB_NUMERIC_COMPARATOR);
    
  /*---Compute-multiplicity-function------------------------------------*/
  Nmax = Ng[Ngroups-1];
  
  indx     = my_calloc(sizeof(*indx),Nmax);
  Ngrps    = my_calloc(sizeof(*Ngrps),Nmax);
  logngrp1 = my_malloc(sizeof(*logngrp1),Nmax);
  logngrp2 = my_malloc(sizeof(*logngrp2),Nmax);

  i=0 ;
  while(Ng[i] < nmin)
    i++;

  assert(i < Ngroups);
  first_index_with_nmin=i;
  for(N=nmin;N<=Nmax;N++) {
    if(i >= Ngroups) {
      fprintf(stderr,"breaking\n");
      break;
    }
    
    if(N==Ng[i]) {
      logngrp1[N-1] = log10((double)(Ngroups-i)*inv_volume);
      while(i < Ngroups && Ng[i]==N) {
	Ngrps[N-1]++ ;
	i++ ;
      }

      if(N==Ng[Ngroups-1]) {
	logngrp2[N-1] = half_log10_inv_volume;
      } else {
	assert(i < Ngroups);	
	logngrp2[N-1] = log10((double)(Ngroups-i)*inv_volume);
      }
    } else if(N<Ng[i]) {
      indx[N-1] = 1 ;
    } else {
      indx[N-1] = 2 ;
    }
  }

  /*---Interpolate-and-extrapolate-------------------------------------*/
  Nlast=first_index_with_nmin;
  for(N=nmin;N<=Nmax;N++) {
    if(indx[N-1]==0) {
      Nlast = N ;
    } else if(indx[N-1]==1) {
      X1 = log10(Nlast) ;
      Y1 = logngrp1[Nlast-1] ;
      Z1 = logngrp2[Nlast-1] ;
      
      j=N ;
      while(j < Nmax && indx[j-1]==1) {
	j++ ;
      }
      assert(j <= Nmax);
      assert(j >= 1);
      
      X2 = log10(j) ;
      Y2 = logngrp1[j-1] ;
      Z2 = logngrp2[j-1] ;
      
      logngrp1[N-1] = Y1 + (log10(N)-X1)*(Y2-Y1)/(X2-X1) ;
      logngrp2[N-1] = Z1 + (log10(N)-X1)*(Z2-Z1)/(X2-X1) ;
    } else {
      X1 = log10(Nlast-1) ;
      Y1 = logngrp1[Nlast-2] ;
      Z1 = logngrp2[Nlast-2] ;
      X2 = log10(Nlast) ;
      Y2 = logngrp1[Nlast-1] ;
      Z2 = logngrp2[Nlast-1] ;
      
      logngrp1[N-1] = Y1 + (log10(N)-X1)*(Y2-Y1)/(X2-X1) ;
      logngrp2[N-1] = Z1 + (log10(N)-X1)*(Z2-Z1)/(X2-X1) ;
    }
  }
  
  /*---Print-output-----------------------------------------------------*/
  
  for(i=nmin-1;i<Nmax;i++) {
    fprintf(stdout,"%3d  %5d  %7.4f %7.4f  %1d\n",i+1,Ngrps[i],logngrp1[i],logngrp2[i],indx[i]);
  }

  free(indx);free(Ngrps);
  free(logngrp1);
  free(logngrp2);
  free(Ng);

  exit(EXIT_SUCCESS);
}
