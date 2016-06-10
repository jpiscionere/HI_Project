#include "redshift_distortions.h"

void pp_zdistortions(float *z, float *vz,float Lbox, const int np, const int zspace,const float zmedian, const int lasdamas_cosmology)
{
  int Nzdc ;
  double *zc,*dc ;
  /*---Redshift-space--------------------*/
  double r,vr;
  /* keeping track of the min and max distortion distance */
  double dd_min, dd_max;
  double dc_max;
	double inv_scale_factor;
  unsigned long nwrap_below=0,nwrap_above=0;
  gsl_interp *interpolation;
  gsl_interp_accel *accelerator;

	fprintf(stderr,"Initializing cosmology with `%s'\n",lasdamas_cosmology==1 ? "LasDamas":"Planck");
	/*Initialize cosmology*/
	init_cosmology(lasdamas_cosmology);
	
  dd_min=9e50;
  dd_max=0;

  if(zmedian < 0) {
#ifdef DEBUG
    fprintf(stderr,"app_pp_zdist> Using particle redshift for distortion distance\n") ;
    fprintf(stderr,"app_pp_zdist> Computing cosmological redshifts\n") ;
#endif
    zc=my_malloc(sizeof(*zc),COSMO_DIST_SIZE);
    dc=my_malloc(sizeof(*dc),COSMO_DIST_SIZE);
    
    Nzdc = set_cosmo_dist(MAX_REDSHIFT_FOR_COSMO_DIST, COSMO_DIST_SIZE, zc, dc,lasdamas_cosmology);
    interpolation = gsl_interp_alloc (gsl_interp_linear,Nzdc);
    accelerator =  gsl_interp_accel_alloc();
    gsl_interp_init(interpolation, dc, zc, Nzdc);//note the reversal -> I want to get redshift for a value of dc (really, r which is a proxy for comoving distance).
    dc_max = dc[Nzdc-1];
  } else {
    /* scale_factor = 1.0/(1+zmedian); */
		inv_scale_factor = (1+zmedian);
    fprintf(stderr,"app_pp_zdist> Using single value for scale factor (single snapshot)\n") ;
  }

  /*---Compute-cosmological-redshifts-for-galaxies----------------------------*/
  double z_local, abs_dd ;
  double hubble, dd ;

  for(int i=0;i<np;i++) {
    r = z[i] + AXIS_OFFSET ; /* offset to start z_min at */
    vr = vz[i];
    if(zmedian < 0) {
      z_local = -1.0;
      assert(r <= dc_max && "co-moving distance is within max distance computed");//make sure that the value of r is inside the populated table
      z_local = gsl_interp_eval(interpolation, dc, zc, r, accelerator); 
      
      if( z_local < 0 ) {
				fprintf(stderr, "add_pp_zdist> ERROR! couldn't set cosmological redshift for particle %d\n", i) ;
				exit(EXIT_FAILURE) ;
      }
      /* scale_factor = 1.0 / (1.0 + z_local) ; */
			inv_scale_factor = (1.0 + z_local);
      /* hubble = HUBBLE * sqrt( OMEGA_M / cube( a ) +  OMEGA_L ); */
      /* dd = vr / a / hubble ; */
    }
    /* hubble = Ho * sqrt( OMEGA_M / cube( scale_factor ) +  OMEGA_L ) ; */
		/* hubble = 100.0 * sqrt( OMEGA_M / cube( scale_factor ) +  OMEGA_L ) ; */
		hubble = 100.0 * sqrt( OMEGA_M *inv_scale_factor*inv_scale_factor*inv_scale_factor +  OMEGA_L ) ;
    /* dd = vr / scale_factor / hubble ; */
		dd = vr * inv_scale_factor / hubble ;
    
    abs_dd = fabs(dd) ;
    
    if(abs_dd < dd_min) {
      dd_min = abs_dd ;
    } else if (abs_dd > dd_max) {
      dd_max = abs_dd ;
    }
    
    z[i] += dd;
    if( z[i] > Lbox ) {
      z[i] -= Lbox ;
      nwrap_above++ ;
    } else if( z[i] < 0 ) {
      z[i] += Lbox ;
      nwrap_below++ ;
    }
  }
  if(zmedian < 0) {
    free(zc);free(dc);
    gsl_interp_free(interpolation);
    gsl_interp_accel_free(accelerator);
    
  } 
#ifdef DEBUG  
  fprintf(stderr,"add_pp_zdist> Distortion distance: \n  min = %g\n  max = %g\n",dd_min, dd_max) ;
  fprintf(stderr,"add_pp_zdist> Numpart wrapped: below = %lu ; above = %lu\n",nwrap_below, nwrap_above) ;
#endif
}

