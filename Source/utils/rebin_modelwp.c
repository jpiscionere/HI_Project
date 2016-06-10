/* PROGRAM rebin_modelwp
   
   --- rebin_modelwp wpfile_obs < wpfile_model > wpfile_model2
   --- re-bins a finely binned model wp(rp) according to an observed wp(rp)

      * wpfile_obs : file containing <rp_min rp_max rp_mean wp(rp) sig_wp> measurements
      * < wpfile_model : file containing <rp xi wp> (output of wp2.c)
      * > wpfile_model2 : file containing <rp wp>
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define MAXBUFSIZE (1000)


int main(int argc, char *argv[])

{
  int i,j ;
/*---arguments-------------------------*/
  FILE *fp=NULL;
/*---wp-files--------------------------*/
  int Nmax,Nbin ;
  double *rp=NULL;
/*---chi2------------------------------*/
  double rpm1,rpm2,wpm1,wpm2,wpmodel ;
  double logrpm1,logrpm2,logwpm1,logwpm2,logwpmodel ;
  char buffer[MAXBUFSIZE];
  
  fp=fopen(argv[1],"r") ;

/*---Read-wp(rp)-observed-file-----------------------------------------*/
  Nmax = 14 ;
  rp=(double *) calloc(Nmax,sizeof(double)) ;

  i=0 ;
  while(1) {
    if(fgets(buffer, MAXBUFSIZE,fp)!=NULL) {
      sscanf(buffer,"%lf ",&rp[i]);
      i++;
    } else {
      break;
    }
  }
  Nbin = i ;
  fclose(fp);
  
/*---Read-wp-model-file-and-rebin--------------------------------------*/
  i=j=0 ;
  rpm1=wpm1=1e-10 ;
  while(1) {
    if(fgets(buffer, MAXBUFSIZE,stdin)!=NULL) {
      sscanf(buffer,"%lf %lf",&rpm2,&wpm2);
      if(rpm2>rp[j]) {
	logrpm1 = log10(rpm1) ;
	logrpm2 = log10(rpm2) ;
	logwpm1 = log10(wpm1) ;
	logwpm2 = log10(wpm2) ;
	
	logwpmodel = logwpm1 + (log10(rp[j])-logrpm1)*(logwpm2-logwpm1)/(logrpm2-logrpm1) ;
	wpmodel = pow(10.,logwpmodel) ;
	
	fprintf(stdout,"%7.4f  %11.5e\n",rp[j],wpmodel);
	j++ ;
      }
      rpm1 = rpm2 ;
      wpm1 = wpm2 ;
      i++ ;
    } else {
      break;
    }
  }
  free(rp);
  return EXIT_SUCCESS;
}


