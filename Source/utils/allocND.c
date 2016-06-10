/* PROGRAM allocND
   ---
   --- Allocation for 2D and 3D arrays of different data types.
   --- usage:
   --- A = alloc2Dint(nx,ny) = alloc2Dfloat(nx,ny) = alloc2Ddouble(nx,ny)
   --- A = alloc3Dint(nx,ny,nz) = alloc3Dfloat(nx,ny,nz) = alloc3Ddouble(nx,ny,nz)
*/

#include <stdio.h>
#include <stdlib.h>
#include "allocND.h"

//-----------------------------------------2D int
int ** alloc2Dint(int nx,int ny)
{
  int i ;
  int **array ;

  array = calloc(nx,sizeof(int *));
  if(array == NULL) {
    fprintf(stderr, "alloc2Dint> out of memory\n");
    exit(EXIT_FAILURE) ;
  }
  for(i = 0; i < nx; i++) {
    array[i] = calloc(ny,sizeof(int));
    if(array[i] == NULL) {
      fprintf(stderr, "alloc2Dint> out of memory\n");
      exit(EXIT_FAILURE) ;
    }
  }
  return array ;
}
//-----------------------------------------2D float
float ** alloc2Dfloat(int nx,int ny)
{
  int i ;
  float **array ;

  array = (float **)calloc(nx,sizeof(float *));
  if(array == NULL) {
    fprintf(stderr, "alloc2Dfloat> Out of memory\n");
    exit(EXIT_FAILURE) ;
  }
  for(i = 0; i < nx; i++) {
    array[i] = (float *)calloc(ny,sizeof(float));
    if(array[i] == NULL) {
      fprintf(stderr, "alloc2Dfloat> out of memory\n");
      exit(EXIT_FAILURE) ;
    }
  }
  return array ;
}
//-----------------------------------------2D double
double ** alloc2Ddouble(int nx,int ny)
{
  int i ;
  double **array ;

  array = (double **)calloc(nx,sizeof(double *));
  if(array == NULL) {
    fprintf(stderr, "alloc2Ddouble> Out of memory\n");
    exit(EXIT_FAILURE) ;
  }
  for(i = 0; i < nx; i++) {
    array[i] = (double *)calloc(ny,sizeof(double));
    if(array[i] == NULL) {
      fprintf(stderr, "alloc2Ddouble> out of memory\n");
      exit(EXIT_FAILURE) ;
    }
  }
  return array ;
}
//-----------------------------------------3D int
int *** alloc3Dint(int nx,int ny,int nz)
{
  int i,j ;
  int ***array ;

  array = calloc(nx,sizeof(int **));
  if(array == NULL) {
    fprintf(stderr, "alloc3Dint> out of memory\n");
    exit(EXIT_FAILURE) ;
  }
  
  for(i = 0; i < nx; i++) {
    array[i] = calloc(ny,sizeof(int *));
    if(array[i] == NULL) {
      fprintf(stderr, "alloc3Dint> out of memory\n");
      exit(EXIT_FAILURE) ;
    }
    for(j = 0; j < ny; j++) {
      array[i][j] = calloc(nz,sizeof(int));
      if(array[i][j] == NULL) {
	fprintf(stderr, "alloc3Dint> out of memory\n");
	exit(EXIT_FAILURE) ;
      }
    }
  }
  return array ;
}

void free2D(void **array,int nx)
{
  for(int i=0;i<nx;i++)
    free(array[i]);

  free(array);
}



//Added by Manodeep
void free3D(void ***array,int nx,int ny)
{
  for(int i=0;i<nx;i++){
    for(int j=0;j<ny;j++){
      free(array[i][j]);
    }
    free(array[i]);
  }
  free(array);
}  



//-----------------------------------------3D float
float *** alloc3Dfloat(int nx,int ny,int nz)
{
  int i,j ;
  float ***array ;

  array = calloc(nx,sizeof(float **));
  if(array == NULL) {
    fprintf(stderr, "alloc3Dfloat> out of memory\n");
    exit(EXIT_FAILURE) ;
  }
  for(i = 0; i < nx; i++) {
    array[i] = calloc(ny,sizeof(float *));
    if(array[i] == NULL) {
      fprintf(stderr, "alloc3Dfloat> out of memory\n");
      exit(EXIT_FAILURE) ;
    }
    for(j = 0; j < ny; j++) {
      array[i][j] = calloc(nz,sizeof(float));
      if(array[i][j] == NULL) {
	fprintf(stderr, "alloc3Dfloat> out of memory\n");
	exit(EXIT_FAILURE) ;
      }
    }
  }
  return array ;
}
//-----------------------------------------3D double
double *** alloc3Ddouble(int nx,int ny,int nz)
{
  int i,j ;
  double ***array ;

  array = calloc(nx,sizeof(double **));
  if(array == NULL) {
    fprintf(stderr, "alloc3Ddouble> out of memory\n");
    exit(EXIT_FAILURE) ;
  }
  for(i = 0; i < nx; i++) {
    array[i] = calloc(ny,sizeof(double *));
    if(array[i] == NULL) {
      fprintf(stderr, "alloc3Ddouble> out of memory\n");
      exit(EXIT_FAILURE) ;
    }
    for(j = 0; j < ny; j++) {
      array[i][j] = calloc(nz,sizeof(double));
      if(array[i][j] == NULL) {
	fprintf(stderr, "alloc3Ddouble> out of memory\n");
	exit(EXIT_FAILURE) ;
      }
    }
  }
  return array ;
}
