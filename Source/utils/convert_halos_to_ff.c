#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ftwrite.h"
#include "utils.h"

//icc -std=c99 convert_halos_to_ff.c utils.c ftwrite.c -xhost -o convert_halos_to_ff

#define MAXLEN   (1000)
int main(int argc, char **argv)
{
  char fname[MAXLEN];
  double *xcen,*ycen,*zcen,*vxcen,*vycen,*vzcen;
  int *pbcflag;
  int *junki;
  long *NumPart;
  double *junkf1,*junkf2;
  int64_t Nhalos;
  int nfiles;
  
  //I need to figure out how to run this code for multiple .halos file (since I need to read in the
  //bgc file to get the header.
  if(argc < 2) { 
    fprintf(stderr,"ERROR: usage `%s'  <list of ascii halo files> \n",argv[0]);
    fprintf(stderr,"argc = %d\n",argc);
    for(int i=1;i<argc;i++)
      fprintf(stderr,"\t\t argv[%d] = %s \n",i,argv[i]);

    exit(EXIT_FAILURE);
  }
  FILE *fp=NULL;
  int nread=0;
  int idat[5];
  float fdat[9];
  float Lbox=1000.0;
  float znow=0.132;
  for(int iarg=1;iarg < argc;iarg++) {
    strncpy(fname,argv[iarg],MAXLEN);
    fprintf(stderr,"Working on file `%s' # %d out of %d files ...",fname,iarg,argc-1);
    Nhalos = getnumlines(fname,'#');
    pbcflag=my_malloc(sizeof(*pbcflag),Nhalos) ;
    NumPart=my_malloc(sizeof(*NumPart),Nhalos) ;
    xcen=my_malloc(sizeof(*xcen),Nhalos) ;
    ycen=my_malloc(sizeof(*ycen),Nhalos) ;
    zcen=my_malloc(sizeof(*zcen),Nhalos) ;
    vxcen=my_malloc(sizeof(*vxcen),Nhalos) ;
    vycen=my_malloc(sizeof(*vycen),Nhalos) ;
    vzcen=my_malloc(sizeof(*vzcen),Nhalos) ;

    //allocate the junk memories
    junki=my_malloc(sizeof(*junki),Nhalos);
    junkf1=my_malloc(sizeof(*junkf1),Nhalos);
    junkf2=my_malloc(sizeof(*junkf2),Nhalos);
    
    fp = my_fopen(fname,"r");

    //read in the data
    for(int i=0;i<Nhalos;i++) {
      nread = fscanf(fp,"%d %ld %lf %lf %lf %lf %lf %lf %lf %lf %d",&junki[i],&NumPart[i],&junkf1[i],&xcen[i],&ycen[i],&zcen[i],&vxcen[i],&vycen[i],&vzcen[i],&junkf2[i],&pbcflag[i]);
      if(nread != 11) {
	fprintf(stderr,"ERROR: Could not parse the string ..exiting\n");
	exit(EXIT_FAILURE);
      }
    }
    fclose(fp);

    //now write out the data in ff format
    my_snprintf(fname,MAXLEN,"%s.ff",argv[iarg]);
    fp = fopen(fname,"r");
    if(fp != NULL) {
      fclose(fp);
      fprintf(stderr,"ERROR: File `%s' already exists..exiting\n",fname);
      exit(EXIT_FAILURE);
    }

    fp=my_fopen(fname,"w");
    idat[0] = (int)Lbox ;
    idat[1] = Nhalos ;
    idat[2]=idat[3]=idat[4]=0;
    fdat[0] = Lbox ;
    fdat[1]=fdat[2]=fdat[3]=fdat[4]=fdat[5]=fdat[6]=fdat[7]=fdat[8]=0.0 ;

    ftwrite(idat,sizeof(int),5,fp);
    ftwrite(fdat,sizeof(float),9,fp);
    ftwrite(&znow,sizeof(znow),1,fp);
    ftwrite(junki,sizeof(*junki),Nhalos,fp);
    ftwrite(NumPart,sizeof(*NumPart),Nhalos,fp);
    ftwrite(junkf1,sizeof(*junkf1),Nhalos,fp);
    ftwrite(xcen,sizeof(*xcen),Nhalos,fp);
    ftwrite(ycen,sizeof(*ycen),Nhalos,fp);
    ftwrite(zcen,sizeof(*zcen),Nhalos,fp);
    ftwrite(vxcen,sizeof(*vxcen),Nhalos,fp);
    ftwrite(vycen,sizeof(*vycen),Nhalos,fp);
    ftwrite(vzcen,sizeof(*vzcen),Nhalos,fp);
    ftwrite(junkf2,sizeof(*junkf2),Nhalos,fp);
    ftwrite(pbcflag,sizeof(*pbcflag),Nhalos,fp);
    fclose(fp);
    
    free(pbcflag);
    free(xcen);free(ycen); free(zcen);
    free(vxcen); free(vycen); free(vzcen);
    free(NumPart);
    free(junki);free(junkf1);free(junkf2);
    fprintf(stderr,"...done\n");
  }

  

  return EXIT_SUCCESS;
}
