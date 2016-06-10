/* PROGRAM wp2
   
   --- wp2 < xifile > wpfile
   --- Computes projected correlation function wp(rp) from the real-space
       correlation function xi(r).  wp(rp)=Int_0^200(dy xi(sqrt(rp^2 + y^2))).
       Extrapolates past maximum scale using a power law.
       Version 2: reads file output by covar3.c rather than xi.c

      * < file containing r, xi (output of covar3.c)
      * > file containing rp, wp
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "/orbital/home/aberlind/Source/Lib/nr.h"
#define sqr(x) ((x)*(x))
#define max(A,B) ((A) > (B) ? (A) : (B))
#define min(A,B) ((A) < (B) ? (A) : (B))
#define mabs(A) ((A) < 0.0 ? -(A) : (A))
#define cnint(x) ((x-floor(x)) < 0.5 ? floor(x) : ceil(x))
#define csign(x) (x < 0.0 ? -1 : 1)
#define PI (3.141592)

int main(int argc, char *argv[])

{
  int i,j ;
/*---xi-file-variables-----------------*/
  int ijunk,Nmax,N ;
  float junk,rmax,*r,*xi ;
  char junkc[15] ;
/*-------------------------------------*/
  float rp,*wp ;
  float y1,y2,y3,y,dy,r1,xi1 ;
  float A,B,C,r0,gamma ;

/*---Read-xi(r)-file---------------------------------------------------*/

  Nmax = 300 ;
  r=(float *) calloc(Nmax,sizeof(float)) ;
  xi=(float *) calloc(Nmax,sizeof(float)) ;

  fscanf(stdin,"%s %s %s %s %s %s %s %s %s %s %s %s %s",&junkc,&junkc,&junkc,&junkc,&junkc,&junkc,&junkc,&junkc,&junkc,&junkc,&junkc,&junkc,&junkc) ;
  i=0 ;
  while(fscanf(stdin,"%f %f %f %f %f %f %f %f %f %f %f %f %f %d",&junk,&junk,&r[i],&xi[i],&junk,&junk,&junk,&junk,&junk,&junk,&junk,&junk,&junk,&ijunk)!=EOF)
    {
      i++ ;
    }
  N = i ;
  rmax = r[N-1] ;

/*---Extrapolate-past-rmax-with-power-law--*/
  gamma = -log10(xi[N-2]/xi[N-3])/log10(r[N-2]/r[N-3]) ;
  C = log10(xi[N-2]) + gamma*log10(r[N-2]) ;
  r0 = pow(10,(C/gamma)) ;
  fprintf(stderr,"wp> Extrapolated power-law: gamma=%f, r0=%f\n",gamma,r0) ;

/*---Compute-wp-for-every-value-of-rp------*/
  wp=(float *) calloc(N,sizeof(float)) ;

  for(i=0;i<N;i++)
    {
      rp = r[i] ;
      
      wp[i] = 0 ;

/*---Integrate-over-y----------------------*/
      y1 = 0 ;
      y2 = sqrt((rmax*rmax - rp*rp)) ;
      y3 = 200 ;
      dy = 0.01 ;
      for(y=y1;y<=y3;y+=dy)
	{
	  r1 = sqrt((rp*rp + y*y)) ;

/*---Interpolate-to-find-xi1=xi(r1)--------*/
	  if(y<y2)
	    {
	      for(j=0;j<N;j++)
		{
		  if(r[j]>r1)
		    {
		      A = log10(xi[j]/xi[j-1]) ;
		      B = log10(r1/r[j-1]) ;
		      C = log10(r[j]/r[j-1]) ;
		      xi1 = pow(10,(log10(xi[j-1]) + A*B/C)) ; 
		      break ;
		    }
		}
	    }
/*---Add-correction-factor-for-y>y2--------*/
	  else
	    {
	      xi1 = pow((r1/r0),(-gamma)) ;
	    }
	  wp[i] += 2*xi1*dy ;
	}

      fprintf(stdout,"%7.4f  %5.3e %5.3e\n",rp,xi[i],wp[i]) ;
    }
  return 0 ;
}


