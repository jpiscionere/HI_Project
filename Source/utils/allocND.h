#ifndef __ALLOCND_H_
#define __ALLOCND_H_

int ** alloc2Dint(int nx,int ny);
float ** alloc2Dfloat(int nx,int ny);
double ** alloc2Ddouble(int nx,int ny);
int *** alloc3Dint(int nx,int ny,int nz);
float *** alloc3Dfloat(int nx,int ny,int nz);
double *** alloc3Ddouble(int nx,int ny,int nz);

void free2D(void **array,int nx);
void free3D(void ***array,int nx,int ny);

#endif
