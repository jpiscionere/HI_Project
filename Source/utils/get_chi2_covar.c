/* PROGRAM get_chi2_covar

   --- get_chi2_covar ngal_obs sig_ngal ngal_model wpfile_obs wpfile_model Nbins diagonal [invert] < covarfile_obs > chi2
   --- Computes chi^2 value from ngal and wp(rp) using covariance matrix.

      * ngal_obs : observed number density of galaxies (h^3/Mpc^3) 
      * sig_ngal : error in observed number density of galaxies (h^3/Mpc^3) 
      * ngal_model : model number density of galaxies (h^3/Mpc^3) 
      * wpfile_obs : file containing <rp_min rp_max rp_mean wp(rp) sig_wp> measurements
      * wpfile_model : file containing <rp wp> (output of rebin_modelwp.c)
      * Nbins : number of bins in wp(rp)
      * diagonal = 0 : use covariance matrix
                 = 1 : use diagonal errors only
      *[invert] = 0 : covariance matrix comes pre-inverted 
                = 1 : invert covariance matrix (default)
      * < covarfile_obs : file containing covariance matrix
      * > chi^2 value
*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#define sqr(x) ((x)*(x))


#define MAXBUFSIZE   (1000)

void Printhelp(void)
{
  fprintf(stderr, "%s", "\n"
  "--- get_chi2_covar ngal_obs sig_ngal ngal_model wpfile_obs wpfile_model Nbins diagonal [invert] < covarfile_obs > chi2\n"
  "--- Computes chi^2 value from ngal and wp(rp) using covariance matrix.\n"
     "* Ngal_obs : observed number density of galaxies (h^3/Mpc^3)\n"
     "* sig_ngal : error in observed number density of galaxies (h^3/Mpc^3)\n"
     "* ngal_model : model number density of galaxies (h^3/Mpc^3)\n" 
     "* wpfile_obs : file containing <rp_min rp_max rp_mean wp(rp) sig_wp> measurements\n"
     "* wpfile_model : file containing <rp wp> (output of rebin_modelwp.c)\n"
     "* Nbins : number of bins in wp(rp)\n"
     "* diagonal = 0 : use covariance matrix\n"
     "           = 1 : use diagonal errors only\n"
     "covarfile_obs : file containing covariance matrix\n"
     "outputfile    : output  chi^2 value\n"
     "*[invert] = 0 : covariance matrix comes pre-inverted\n"
     "          = 1 : invert covariance matrix (default)\n"
	  );
}

int main(int argc, char *argv[])
{
  int     i,j,l,input_format,diagonal_errors;
  FILE    *fp1=NULL,*fp2=NULL;
  char    *data_file=NULL,*theory_file=NULL,*covar_file=NULL,*output_file=NULL;
  void    inv (double**,int);
  double  chisquared_function (int , double**, double*, double*, double*, double);

  int     N_file;
  double  *datawp=NULL;
  double  chisquared;  
  //double  fjunk;
  int     bins2=14;
  double  *theory=NULL,*error=NULL;
  double  **covariance_matrix=NULL;
  double  num_den,err_num_den,theory_num_den;
  char buffer[MAXBUFSIZE];
  const char argnames[][100] = {"Ngal (obs)","Sig Ngal (obs)","Ngal (model)","wp file (obs)","wp file(model)","Nbins","Diagonal (0-cov. matrix,1=diag. errors)",
				"Covar file (obs)","Output file","[invert Cov. matrix ()]"};
  int nargs=sizeof(argnames)/(sizeof(char)*100);
//---Read-arguments-----------------------//
  if(argc<10)  {
      Printhelp() ;
      return EXIT_FAILURE ;
    }
  sscanf(argv[1],"%lf",&num_den);
  sscanf(argv[2],"%lf",&err_num_den);
  sscanf(argv[3],"%lf",&theory_num_den);
  data_file = argv[4];
  theory_file = argv[5];
  sscanf(argv[6],"%d",&bins2);
  sscanf(argv[7],"%d",&diagonal_errors);
  covar_file=argv[8];
  output_file=argv[9];

  input_format = 1 ;
  if(argc>10) sscanf(argv[10],"%d",&input_format);


  fprintf(stderr,"Running `%s' with the parameters \n",argv[0]);
  fprintf(stderr,"\n\t\t-----------------------------------------------------\n");
  for(i=1;i<argc;i++) {
    if(i <= nargs) {
      fprintf(stderr,"\t\t %-45s = %s \n",argnames[i-1],argv[i]);
    }  else {
      fprintf(stderr,"\t\t <> = `%s' \n",argv[i]);
    }
  }
  fprintf(stderr,"\t\t-----------------------------------------------------\n\n");
  
    
  
//---Read-wp-data-------------------------//
  error=(double *) calloc(bins2,sizeof(double)) ;
  datawp=(double *) calloc(bins2,sizeof(double)) ;
  assert(error != NULL);
  assert(datawp != NULL);
  
  fp1=fopen(data_file,"r") ;
  assert(fp1 != NULL );
  
  l=0;
  while(1){
    if(fgets(buffer, MAXBUFSIZE,fp1)!=NULL) {
      //sscanf(buffer,"%lf %lf %lf %lf %lf",&fjunk,&fjunk,&fjunk,&datawp[l],&error[l]);
      sscanf(buffer,"%*lf %*lf %*lf %lf %lf",&datawp[l],&error[l]);//suppress assignment
      l++;
    } else {
      break;
    }
  }

  N_file = l ;
  
  if(N_file != bins2) {
    fprintf(stderr,"covariance > BIN NUMBER MISMATCH!! You DIDN'T filter the DATA FILE `%s' (N_file=%d, bins2=%d)!\n",data_file,N_file,bins2);
    fprintf(stderr,"covariance > BAILING OUT!\n");
    fprintf(stderr,"N_file=%d\n",N_file);
    return EXIT_FAILURE;		
  }
  
  fclose(fp1) ;
  
//---Read-wp-theory-----------------------//

  theory=(double *) calloc(bins2,sizeof(double));
  fp2=fopen(theory_file,"r") ;
  assert(fp2 != NULL );
  assert(theory != NULL);
  
  l=0;
  while(1){
    if(fgets(buffer, MAXBUFSIZE,fp2)!=NULL) {
      //sscanf(buffer,"%lf %lf ",&fjunk,&theory[l]);
      sscanf(buffer,"%lf ",&theory[l]);
      l++;
    } else {
      break;
    }
  }
  
  N_file = l ;
  if(N_file != bins2) {
    fprintf(stderr,"covariance > BIN NUMBER MISMATCH!! You DIDN'T filter the THEORY FILE `%s' ( N_file = %d, bins2 = %d)!\n",theory_file,N_file,bins2);
    fprintf(stderr,"covariance > BAILING OUT!\n");
    return EXIT_FAILURE;                
  }
  fclose(fp2) ;	
  
//---Read-covar-data----------------------//
  covariance_matrix = (double **) malloc(bins2 * sizeof(double *));
  assert(covariance_matrix != NULL);
  for (i = 0; i < bins2; i++) {
    covariance_matrix[i] = (double *) malloc(bins2 * sizeof(double));
    assert(covariance_matrix[i] != NULL);
  }

  fp1 = fopen(covar_file,"r");
  assert(fp1 != NULL);

  for(i=0;i<bins2;i++) {
    /* if(fgets(buffer,MAXBUFSIZE,fp1)==NULL) { */
    /*   fprintf(stderr,"get_chi2_covar> ERROR - could not read covariance - matrix..exiting\n"); */
    /*   exit(EXIT_FAILURE); */
    /* } */
    /* fprintf(stderr,"buffer = %s \n",buffer); */
    for(j=0;j<bins2;j++) {
      fscanf(fp1,"%lf ",&covariance_matrix[i][j]);
    }


    /* for(j=0;j<bins2;j++) { */
    /*   fprintf(stderr,"%lf ",covariance_matrix[i][j]); */
    /* } */
    /* fprintf(stderr,"\n"); */
  }
  fclose(fp1);
  
//---Invert-matrix------------------------//

  if(input_format == 1)
    inv(covariance_matrix,bins2);	
  
//---Compute-chi^2------------------------//

  chisquared=chisquared_function(bins2,covariance_matrix, theory, datawp, error, diagonal_errors);
  {
    double intermediate_chisqr= sqr((num_den-theory_num_den)/(err_num_den));
#ifdef PRINT_INDIVIDUAL_CHISQR
    fprintf(stderr,"%3d  %20.10lf \n",bins2,intermediate_chisqr);
#endif      
    chisquared += intermediate_chisqr;
  }

  
  fprintf(stderr,"get_chi2_covar> chi^2 = %lf \n",chisquared);

  //check that the output file is not about to over-write something else (i.e., covar file)
  fp2=fopen(output_file,"r");
  assert(fp2==NULL);
  /* fclose(fp2); */ ///not reqd. since the file should not exist

  //now open the output file for realsies
  fp2 = fopen(output_file,"w");
  assert(fp2 != NULL);
  fprintf(fp2,"%lf \n",chisquared);
  fclose(fp2);

  free(error);
  free(datawp);
  free(theory);
  for (i = 0; i < bins2; i++) {
    free(covariance_matrix[i]);
  }
  free(covariance_matrix);
  
  return EXIT_SUCCESS;
}

void inv(double **a_matrix, int dimension)
{  
  int i,j;
  
  gsl_matrix * C = gsl_matrix_calloc ( dimension, dimension );
  gsl_matrix * C_i = gsl_matrix_calloc ( dimension, dimension );
  
  for ( i=0;i<dimension;i++ ) {
    for ( j=0;j<dimension;j++ ) {
      gsl_matrix_set ( C, i, j, a_matrix[i][j] );
    }
  }
  int s;
  
  gsl_permutation *p = gsl_permutation_alloc ( dimension );
  
  gsl_linalg_LU_decomp ( C, p, &s );
  
  gsl_linalg_LU_invert ( C, p, C_i );
  
  for ( i=0;i<dimension;i++ ) {
    for ( j=0;j<dimension;j++ ) {
      a_matrix[i][j] = gsl_matrix_get( C_i, i, j ) ;
      
    }
  }  
  gsl_permutation_free( p );
  gsl_matrix_free( C );
  gsl_matrix_free( C_i );
}

double chisquared_function (int dimen, double **covar, double *theory, double *data, double *error, double diagonal_errors)
{
  double *intermediate_step=NULL,chisquared=0.0,intermediate_chisqr=0.0;
  int i,j;

  intermediate_step = (double * ) calloc(dimen,sizeof(double));
  assert(intermediate_step != NULL);
  for(i=0;i<dimen;i++) {
    intermediate_step[i] = (data[i]-theory[i])/(error[i]);    
  }
  
  for(i=0;i<dimen;i++){
    if(diagonal_errors==0) {
      intermediate_chisqr = 0.0;
      for(j=0;j<dimen;j++) {
	intermediate_chisqr+=intermediate_step[i]*covar[i][j]*intermediate_step[j];//use covariance matrix
      }
    } else {
      intermediate_chisqr = intermediate_step[i]*intermediate_step[i]; //only use diagonal errors
    }
    chisquared +=intermediate_chisqr;

#ifdef PRINT_INDIVIDUAL_CHISQR 
      fprintf(stderr,"%3d  %20.10lf \n",i,intermediate_chisqr);
#endif      
      
  }
  
  free(intermediate_step);
  return chisquared;  
} 
