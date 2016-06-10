#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "ftread.c"

int main(int argc, char **argv)
{
  if(argc < 2) {
    fprintf(stderr,"ERROR: Usage `%s'  <fast-food files>\n",argv[0]);
    exit(EXIT_FAILURE);
  }

  for(int iarg=1;iarg<argc;iarg++) {
    FILE *fp=fopen(argv[iarg],"r");
    assert(fp != NULL);
    int idat[5];
    float fdat[9];
    my_ftread(idat,sizeof(*idat),5,fp);
    my_ftread(fdat,sizeof(*fdat),9,fp);
    fclose(fp);
    fprintf(stderr,"Running `%s' on the file `%s'\n\n",argv[0],argv[iarg]);
    
    fprintf(stderr," BoxSize = %12.2lf\n Ngal    = %12d\n ngal    = %12.6lf \n\n",fdat[0],idat[1],idat[1]/(fdat[0]*fdat[0]*fdat[0]));
    if(argc == 2) {
      fprintf(stderr,"***** HOD ***** \n");
      fprintf(stderr," Ncenopt = %2d \n Nsatopt = %2d \n PNNopt  = %2d\n",idat[2],idat[3],idat[4]);
      fprintf(stderr,"*************** \n\n");
      
      fprintf(stderr,"++++ Parameter Values ++++ \n");
      fprintf(stderr," logMmin = %12.6lf \n siglogM = %12.6lf \n logM0   = %12.6lf \n logM1   = %12.6lf \n alpha   = %12.6lf \n fsat    = %12.6f \n fgal    = %12.6f \n gamma   = %12.6f \n",
	      fdat[1],fdat[2],fdat[3],fdat[4],fdat[5],fdat[6],fdat[7],fdat[8]);
      fprintf(stderr,"++++++++++++++++++++++++++ \n\n");
    }

  }
    
  exit(EXIT_SUCCESS);
}
