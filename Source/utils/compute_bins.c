#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//Create the bins for input to wp_covar

int main(int argc,char **argv)
{
  double xmin,xmax,binsize;
  int Nbins,logbins;
  if(argc != 5) {
    fprintf(stderr,"ERROR: usage - `%s' <xmin> <xmax>  <Nbins>  <log=1,linear=otherwise>\n",argv[0]);
    exit(EXIT_FAILURE);
  }

  xmin  = atof(argv[1]);
  xmax  = atof(argv[2]);
  Nbins = atoi(argv[3]);
  logbins = atoi(argv[4]);

  //error-check parameters
  if(xmin == xmax) {
    fprintf(stderr,"ERROR: xmin can not be equal to xmax ..exiting\n");
    exit(EXIT_FAILURE);
  }

  if(Nbins <= 0) {
    fprintf(stderr,"ERROR: Nbins must be greater than 0 ..exiting\n");
    exit(EXIT_FAILURE);
  }

  
  fprintf(stderr,"Running `%s' with the parameters \n",argv[0]);
  fprintf(stderr,"\t\t xmin    = %lf\n",xmin);
  fprintf(stderr,"\t\t xmax    = %lf\n",xmax);
  fprintf(stderr,"\t\t Nbins   = %d \n",Nbins);
  fprintf(stderr,"\t\t logbins = %d \n",logbins);
  
  if(logbins == 1) {
    //create logbins
    //check that xmin and xmax are > 0
    if(xmin < 0.0 || xmax < 0.0) {
      fprintf(stderr,"ERROR: Can not make logarithmic bins with negative xrange <xmin,xmax> = <%lf,%lf>..exiting\n",xmin,xmax);
      exit(EXIT_FAILURE);
    }

    binsize = log(xmax/xmin)/Nbins;
    for(int i=0;i<Nbins;i++)
      fprintf(stdout,"%lf %lf \n",exp(i*binsize+log(xmin)),exp((i+1)*binsize+log(xmin)));
  } else {
    binsize = (xmax-xmin)/Nbins;
    for(int i=0;i<Nbins;i++)
      fprintf(stdout,"%lf %lf \n",i*binsize+xmin,(i+1)*binsize+xmin);
  }
  
  return EXIT_SUCCESS;
}
