/* PROGRAM get_joint_chi2_covar*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#include "utils.h"

/* #define sqr(x) ((x)*(x)) */

#define MAXLEN       (1000)
#define MAXBUFSIZE   (10000)

/* void Printhelp(void) */
/* { */
/*   fprintf(stderr, "%s", "\n" */
/*   "--- get_chi2_covar ngal_obs sig_ngal ngal_model wpfile_obs wpfile_model Nbins diagonal [invert] < covarfile_obs > chi2\n" */
/*   "--- Computes chi^2 value from ngal and wp(rp) using covariance matrix.\n" */
/*      "* Ngal_obs : observed number density of galaxies (h^3/Mpc^3)\n" */
/*      "* sig_ngal : error in observed number density of galaxies (h^3/Mpc^3)\n" */
/*      "* ngal_model : model number density of galaxies (h^3/Mpc^3)\n"  */
/*      "* wpfile_obs : file containing <rp_min rp_max rp_mean wp(rp) sig_wp> measurements\n" */
/*      "* wpfile_model : file containing <rp wp> (output of rebin_modelwp.c)\n" */
/*      "* Nbins : number of bins in wp(rp)\n" */
/*      "* diagonal = 0 : use covariance matrix\n" */
/*      "           = 1 : use diagonal errors only\n" */
/*      "covarfile_obs : file containing covariance matrix\n" */
/*      "outputfile    : output  chi^2 value\n" */
/*      "*[invert] = 0 : covariance matrix comes pre-inverted\n" */
/*      "          = 1 : invert covariance matrix (default)\n" */
/* 	  ); */
/* } */

int main(int argc, char *argv[])
{
  int     i,j,l,diagonal_errors=0;
  FILE    *fp1=NULL,*fp2=NULL;
  //char    *data_file=NULL,*theory_file=NULL,*covar_file=NULL,*output_file=NULL;
  void    inv (double**,int);
  //double  chisquared_function (int , double**, double*, double*, double*, double);
  double chisquared_function (int dimen, double **covar, double *chi, int diagonal_errors);
  double chisqr_pca(const int totnbins, double **correlation_matrix,const double threshold,double *chi);
  int     N_file;
  double  *data=NULL;
  double  chisquared;  
  //double  fjunk;
  /* int     bins2=14; */
  double  *theory=NULL,*error=NULL;
  double  **correlation_matrix=NULL,**effective_correlation_matrix=NULL;
  double  num_den,err_num_den,theory_num_den;
  char buffer[MAXBUFSIZE];
  int fit_wp,fit_gmf,wp_Nbins,gmf_Nbins;
  char sdss_ngal_file[MAXLEN],sdss_wp_file[MAXLEN],sdss_gmf_file[MAXLEN];
  char wp_file[MAXLEN],gmf_file[MAXLEN];
  char covar_file[MAXLEN],output_file[MAXLEN];
  double Ngal_theory;
  //my_snprintf(execstring,MAXLEN,"%s/get_joint_chi2_covar %d %d  %s %d %d %s %s %s %s %s %s",bindir,fit_wp,fit_gmf,sdss_ngal_file,wp_Nbins,gmf_Nbins,sdss_wp_file,wp_file,sdss_gmf_file,gmf_file_base,sdss_joint_cov_file,chi2_file);
  
  const char argnames[][100] = {"wp fit","gmf fit","SDSS Ngal file","Nbins (wp)","Nbins (gmf)","wp file (SDSS)","wp file (Model)","Ngal (Model)","gmf file(SDSS)","gmf file(Model)",
				"Covar file (Joint)","Output file","Use PCA","Threshold (for trimming modes in PCA)"};
  int nargs=sizeof(argnames)/(sizeof(char)*100);
  const int optional_input=0;
  int required_nargs = nargs-optional_input +1;//the "+1" is for argv[0]

  //---Read-arguments-----------------------//
  if(argc<required_nargs)  {
    //Printhelp() ;
    fprintf(stderr,"ERROR: usage: `%s' \n",argv[0]);
    if(argc > 1) {
      fprintf(stderr,"\n\t\t Found parameters:\n");
      fprintf(stderr,"\t\t --------------------------------------------\n");
      for(i=1;i<argc;i++) {
	fprintf(stderr,"\t\t %-30s = %s \n",argnames[i-1],argv[i]);
      }
      fprintf(stderr,"\t\t --------------------------------------------\n");
    }
    
    fprintf(stderr,"\n\t\t Missing parameters:\n");
    fprintf(stderr,"\t\t *********************************************\n");
    for(i=argc;i<=nargs;i++) {
      fprintf(stderr,"\t\t %-30s = `?'\n",argnames[i-1]);
    }
    fprintf(stderr,"\t\t *********************************************\n");
    return EXIT_FAILURE ;
  }
  fit_wp=atoi(argv[1]);
  fit_gmf=atoi(argv[2]);
  my_snprintf(sdss_ngal_file,MAXLEN,"%s",argv[3]);
  wp_Nbins=atoi(argv[4]);
  gmf_Nbins=atoi(argv[5]);

  my_snprintf(sdss_wp_file,MAXLEN,"%s",argv[6]);
  my_snprintf(wp_file,MAXLEN,"%s",argv[7]);
  Ngal_theory = atof(argv[8]);
  my_snprintf(sdss_gmf_file,MAXLEN,"%s",argv[9]);
  my_snprintf(gmf_file,MAXLEN,"%s",argv[10]);
  my_snprintf(covar_file,MAXLEN,"%s",argv[11]);
  my_snprintf(output_file,MAXLEN,"%s",argv[12]);
  int use_pca = atoi(argv[13]);
  double threshold = atof(argv[14]);
  
  /* input_format = 1 ; */
  /* if(argc>=required_nargs) sscanf(argv[required_nargs-1],"%d",&input_format); */


  fprintf(stderr,"Running `%s' with the parameters \n",argv[0]);
  fprintf(stderr,"\n\t\t-----------------------------------------------------\n");
  for(i=1;i<argc;i++) {
    if(i <= nargs) {
      fprintf(stderr,"\t\t %-30s = %s \n",argnames[i-1],argv[i]);
    }  else {
      fprintf(stderr,"\t\t <> = `%s' \n",argv[i]);
    }
  }
  fprintf(stderr,"\t\t-----------------------------------------------------\n\n");
  
  int covar_Nbins = 1 + wp_Nbins + gmf_Nbins ; //The first bin is for Ngal, next wp_Nbins is for wp; next gmf_Nbins is for gmf
  int TotNbins=1;//for Ngal.   
  int lstart;
  int nread,nitems;

  if(fit_wp==1)
    TotNbins += wp_Nbins;

  if(fit_gmf==1)
    TotNbins += gmf_Nbins;
  
  
  //---Read-wp-data-------------------------//
  error  = my_malloc(sizeof(*error),TotNbins);
  data   = my_malloc(sizeof(*data),TotNbins);
  theory = my_malloc(sizeof(*theory),TotNbins);
  
  //fill in Ngal + sig_Ngal info
  fp1=my_fopen(sdss_ngal_file,"r") ;
  fscanf(fp1,"%lf %lf",&data[0],&error[0]);
  fclose(fp1);

  //Okay - now fill in the theory measurements
  theory[0] = Ngal_theory;
  
  if(fit_wp==1) {
    lstart=1;
    fp1=my_fopen(sdss_wp_file,"r") ;
    
    l=lstart;
    nitems=2;
    while(fgets(buffer, MAXBUFSIZE,fp1)!=NULL){
      //sscanf(buffer,"%lf %lf %lf %lf %lf",&fjunk,&fjunk,&fjunk,&datawp[l],&error[l]);
      nread=sscanf(buffer,"%*lf %*lf %*lf %lf %lf",&data[l],&error[l]);//suppress assignment
      assert(nread==nitems);
      l++;
    }
    N_file = l-lstart ;
    
    if(N_file != wp_Nbins) {
      fprintf(stderr,"covariance > BIN NUMBER MISMATCH!! You DIDN'T filter the DATA FILE `%s' (N_file=%d, wp_nbins=%d)!\n",sdss_wp_file,N_file,wp_Nbins);
      fprintf(stderr,"covariance > BAILING OUT!\n");
      return EXIT_FAILURE;		
    }
    fclose(fp1) ;

    fp2=my_fopen(wp_file,"r") ;
    l=lstart;
    nitems=1;
    while(fgets(buffer, MAXBUFSIZE,fp2)!=NULL) {
      nread=sscanf(buffer,"%lf ",&theory[l]);
      assert(nitems==nread);
      l++;
    }
    N_file = l-lstart ;
    if(N_file != wp_Nbins) {
      fprintf(stderr,"covariance > BIN NUMBER MISMATCH!! You DIDN'T filter the MODEL FILE `%s' (N_file=%d, wp_nbins=%d)!\n",wp_file,N_file,wp_Nbins);
      fprintf(stderr,"covariance > BAILING OUT!\n");
      return EXIT_FAILURE;                
    }
    fclose(fp2) ;	
  }


  if(fit_gmf==1) {
    fp1=my_fopen(sdss_gmf_file,"r") ;

    if(fit_wp==1) {
      lstart=wp_Nbins+1;
    } else {
      lstart=1;
    }

    l=lstart;
    nitems=2;
    while(fgets(buffer, MAXBUFSIZE,fp1)!=NULL){
      //sscanf(buffer,"%lf %lf %lf %lf %lf",&fjunk,&fjunk,&fjunk,&datawp[l],&error[l]);
      nread=sscanf(buffer,"%*d %*d %lf %lf",&data[l],&error[l]);//suppress assignment
      assert(nread==nitems);
      l++;
    }
    N_file = l-lstart ;
    
    if(N_file != gmf_Nbins) {
      fprintf(stderr,"covariance > BIN NUMBER MISMATCH!! You DIDN'T filter the DATA FILE `%s' (N_file=%d, gmf_nbins=%d)!\n",sdss_gmf_file,N_file,gmf_Nbins);
      fprintf(stderr,"covariance > BAILING OUT!\n");
      return EXIT_FAILURE;		
    }
    fclose(fp1) ;

    //Okay - now fill in the theory measurements
    fp2=my_fopen(gmf_file,"r") ;
    
    l=lstart;
    nitems=1;
    while(fgets(buffer, MAXBUFSIZE,fp2)!=NULL) {
      nread=sscanf(buffer,"%lf ",&theory[l]);
      assert(nread==nitems);
      l++;
    }
    N_file = l-lstart ;
    if(N_file != gmf_Nbins) {
      fprintf(stderr,"covariance > BIN NUMBER MISMATCH!! You DIDN'T filter the MODEL FILE `%s' (N_file=%d, wp_nbins=%d)!\n",gmf_file,N_file,gmf_Nbins);
      fprintf(stderr,"covariance > BAILING OUT!\n");
      return EXIT_FAILURE;                
    }
    fclose(fp2) ;	
  }

  /* for(i=0;i<TotNbins;i++) */
  /*   fprintf(stderr,"Theory = %g data = %g error = %g \n",theory[i],data[i],error[i]); */
  
  //---Read-entire-covar-data----------------------//
  correlation_matrix = (double **) my_malloc(sizeof(double *),covar_Nbins);
  for (i = 0; i < covar_Nbins; i++) {
    correlation_matrix[i] = (double *) my_malloc(sizeof(double),covar_Nbins);
  }

  fp1 = my_fopen(covar_file,"r");
  nitems=1;
  for(i=0;i<covar_Nbins;i++) {
    for(j=0;j<covar_Nbins;j++) {
      nread=fscanf(fp1,"%lf ",&correlation_matrix[i][j]);
      assert(nread==nitems);
    }
  }
  fclose(fp1);

  //Now assign the covariance matrix to the appropriate place
  effective_correlation_matrix = (double **) my_malloc(sizeof(double *),TotNbins);
  for (i = 0; i < TotNbins; i++) {
    effective_correlation_matrix[i] = (double *) my_malloc(sizeof(double),TotNbins);
  }

  if(fit_wp==1 && fit_gmf==1) {
    assert(covar_Nbins == TotNbins);
    for(i=0;i<TotNbins;i++) {
      for(j=0;j<TotNbins;j++) {
	effective_correlation_matrix[i][j]=correlation_matrix[i][j];
      }
    }
  }

  if(fit_wp==1 && fit_gmf==0) {
    for(i=0;i<TotNbins;i++) {
      for(j=0;j<TotNbins;j++) {
	effective_correlation_matrix[i][j]=correlation_matrix[i][j];
      }
    }
  }

  if(fit_wp==0 && fit_gmf==1) {
    int xindex,yindex;
    for(i=0;i<TotNbins;i++) {
      for(j=0;j<TotNbins;j++) {
	xindex = i + wp_Nbins ;
	yindex = j + wp_Nbins ;
	if(i==0 && j==0) {
	  xindex=i;
	  yindex=j;
	}

	if(i==0 && j !=0) {
	  xindex=i;
	  yindex=j+wp_Nbins;
	}

	if(j==0 && i != 0) {
	  xindex = i+wp_Nbins;//index: 0->Ngal;(1-wp_Nbins)->wp;(wp_Nbins+1:wp_Nbins+gmf_Nbins+1) -> gmf
	  yindex=j;
	}

	assert(xindex < covar_Nbins);
	assert(yindex < covar_Nbins);

	effective_correlation_matrix[i][j]=correlation_matrix[xindex][yindex];
      }
    }
  }


  /* for(i=0;i<TotNbins;i++) { */
  /*   for(j=0;j<TotNbins;j++) { */
  /*     fprintf(stderr,"%8.3lf ",effective_correlation_matrix[i][j]); */
  /*   } */
  /*   fprintf(stderr,"\n"); */
  /* } */

  
  
  /* //---Invert-matrix------------------------// */
  /* if(input_format == 1) { */
  /*   inv(effective_correlation_matrix,TotNbins); */
  /* } */
  
  double *intermediate_step=NULL;
  intermediate_step = my_malloc(sizeof(*intermediate_step), TotNbins);
  for(i=0;i<TotNbins;i++) {
    intermediate_step[i] = (data[i]-theory[i])/(error[i]);    
  }

  //---Compute-chi^2------------------------//
  chisquared=0.0;
  if(use_pca == 0) {
    chisquared = chisquared_function(TotNbins, effective_correlation_matrix, intermediate_step, diagonal_errors);
  } else {
    chisquared = chisqr_pca(TotNbins,effective_correlation_matrix, threshold, intermediate_step);
  }
  free(intermediate_step);
  fprintf(stderr,"get_chi2_covar> chi^2 = %lf \n",chisquared);

  //check that the output file is not about to over-write something else (i.e., covar file)
  fp2=fopen(output_file,"r");
  if(fp2==NULL) {
  /* fclose(fp2); */ ///not reqd. since the file should not exist

  //now open the output file for realsies
    fp2 = fopen(output_file,"w");
    assert(fp2 != NULL);
    fprintf(fp2,"%lf \n",chisquared);
    fclose(fp2);
  } else {
    fclose(fp2);
    fprintf(stderr,"Not writing out chi2 file since it already exists\n");
  }
  free(error);
  free(data);
  free(theory);
  for (i = 0; i < TotNbins; i++) {
    free(effective_correlation_matrix[i]);
  }
  free(effective_correlation_matrix);

  for (i = 0; i < covar_Nbins; i++) {
    free(correlation_matrix[i]);
  }
  free(correlation_matrix);

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


double chisqr_pca(const int totnbins, double **correlation_matrix,const double threshold,double *chi)
{
  gsl_matrix *m = gsl_matrix_alloc(totnbins,totnbins);
  for(int i=0;i<totnbins;i++) {
    for(int j=0;j<totnbins;j++) {
      gsl_matrix_set(m, i, j, correlation_matrix[i][j]);
    }
  }

  gsl_vector *eval = gsl_vector_alloc (totnbins);
  gsl_matrix *evec = gsl_matrix_alloc (totnbins, totnbins);
  gsl_eigen_symmv_workspace * w  = gsl_eigen_symmv_alloc(totnbins);
  gsl_eigen_symmv(m, eval, evec, w);
  gsl_eigen_symmv_free (w);
  gsl_matrix_free(m);
  gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_VAL_DESC );

  int L=0;
  for (int i = 0; i < totnbins; i++) {
    double eval_i = gsl_vector_get(eval, i);
    gsl_vector_view evec_i = gsl_matrix_column(evec, i);

    if(eval_i >= threshold) {
      /* fprintf(stdout,"PASS: i = %d eigenvalue = %10.6lf\n", i+1,eval_i); */
      L++;
    }
  }
  assert(L > 0);
  gsl_vector_view vector_chi = gsl_vector_view_array(chi, totnbins);
  gsl_vector *temp_vector = gsl_vector_alloc(totnbins);
  gsl_vector *projected  = gsl_vector_alloc(L);
  gsl_matrix *pc_space   = gsl_matrix_alloc(totnbins, L);
  for (int i=0;i<L; i++){
    gsl_matrix_get_col(temp_vector, evec, i);
    gsl_matrix_set_col(pc_space, i, temp_vector);
  }
  gsl_blas_dgemv(CblasTrans , 1.0, pc_space, &vector_chi.vector, 0.0, projected);
  gsl_vector_free(temp_vector);
  gsl_matrix_free(pc_space);

  double chisqr_val=0.0;
  double eval_i,projected_i;
  for(int i=0;i < L;i++) {
    projected_i = gsl_vector_get(projected, i);
    eval_i = gsl_vector_get(eval, i);
    assert(eval_i >= threshold);
    chisqr_val += projected_i*projected_i/eval_i;
    /* fprintf(stderr,"PROJECTED: i = %d , %10.6lf chi2 = %10.6lf \n",i+1,projected_i,projected_i*projected_i/eval_i); */
  }
  fprintf(stderr,"chisqr_pca> DONE: chisqr = %10.6lf L = %d\n",chisqr_val,L);
  gsl_matrix_free(evec);
  gsl_vector_free(eval);
  gsl_vector_free(projected);

  return chisqr_val;
}



double chisquared_function (int dimen, double **covar, double *intermediate_step, int diagonal_errors)
{
  double chisquared=0.0,intermediate_chisqr=0.0;
  int i,j;

  //Correlation matris is Always supplied - inverted here. 
  inv(covar,dimen);

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
  

  return chisquared;  
} 
