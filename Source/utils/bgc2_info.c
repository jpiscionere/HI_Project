#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <assert.h>
#include <inttypes.h>


#include "bgc2_read_utils.h"

int main( int argc, char **argv )
{
  OUTPUT_HEADER hdr;
  FILE *fp ;
  int nmin_for_group=20;

  if( argc < 2 ) {
    fprintf( stderr, "Usage: %s BGC_FILE1 BGC_FILE2 ... \n", argv[0] );
    exit( EXIT_FAILURE );
  }
  
  // First read and write header
  for(int i=1;i<argc;i++) {
	fp = fopen( argv[i], "r" );
	assert( fp != 0 );
	bgc_read_header( fp, &hdr );
	fprintf(stderr,"********* File = `%s' ************* \n",argv[i]);
	fprintf(stderr,"Total Groups = %"PRId64"\n",hdr.ngroups_total);
	fprintf(stderr,"Particle mass = %lf \n",hdr.part_mass);
	fprintf(stderr,"Min. particles per group = %"PRId64" log10 min. group mass = %lf \n",hdr.min_group_part,log10(hdr.min_group_part*hdr.part_mass));
	fprintf(stderr,"log10 Min mass for %d particles = %lf\n",nmin_for_group,log10(hdr.part_mass*nmin_for_group));
	fprintf(stderr,"*************************************\n\n");
	fclose(fp);
  }
  exit(EXIT_SUCCESS);
}
