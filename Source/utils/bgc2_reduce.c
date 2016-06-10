#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <assert.h>
#include <inttypes.h>

/* this utility reads a set of BGC2 files with particle info and outputs a single BGC2 file with only header and group info
 * */

/* #include "bgc2_read_utils.c" */
/* #include "bgc2_write_utils.c" */

#include "bgc2_read_utils.h"
#include "bgc2_write_utils.h"

static int64_t zero_sat_count=0;
void group_readwrite( char *filename )
{
    FILE *fp;
    OUTPUT_HEADER hdr;
    void *gdata;
    int64_t nhosts=0,nsats=0;

    fprintf(stderr,"Reading '%s' ...",filename);
    fp = fopen( filename, "r" );
    assert( fp != 0 );

    bgc_read_header( fp, &hdr );
    gdata = bgc_read_grouplist( fp, hdr );
    bgc_write_grouplist_raw( stdout, hdr.ngroups, hdr.format_group_data, gdata ) ;
    for(int64_t i=0;i<hdr.ngroups;i++) {
      GROUP_DATA_RMPVMAX gd;
      gd = ((GROUP_DATA_RMPVMAX  *)gdata)[i];
      if(gd.parent_id >=0 ) {
	nsats++;
      } else {
	nhosts++;
      }
    }
    if(nsats==0)
      zero_sat_count++;
    
    fprintf(stderr," ngroups = %"PRId64" nhosts = %"PRId64" nsats = %"PRId64" sat. frac = %e\n",
	    hdr.ngroups,nhosts,nsats,nsats/(double) hdr.ngroups);
    fclose( fp );
    free(gdata);
}

int main( int argc, char **argv )
{
  OUTPUT_HEADER hdr;
  FILE *fp ;
  int64_t ngroups_total;
  int64_t npart_total;
  int64_t max_npart=0;
  
  if( argc < 2 ) {
    fprintf( stderr, "Usage: %s BGC_FILE1 BGC_FILE2 ... > BGC_output\n", argv[0] );
    exit( EXIT_FAILURE );
  }
  
  // First read and write header
  fp = fopen( argv[1], "r" );
  assert( fp != NULL );
  bgc_read_header( fp, &hdr );
  fclose(fp);
  fprintf(stderr,"Reading file `%s' ngroups = %"PRId64"\n",argv[1],hdr.ngroups);
  if(hdr.num_files != argc-1) {
    fprintf(stderr,"ERROR: %ld files expected, only %d found!\n",(ssize_t)hdr.num_files,(argc-1)) ;
    exit(EXIT_FAILURE);
  }

  if(hdr.ngroups_total == 0) {
    ngroups_total = 0;
    for(int i=1; i<argc; i++ ) {
      fp = fopen(argv[i],"r");
      assert( fp != NULL);
      bgc_read_header( fp, &hdr );
      ngroups_total += hdr.ngroups;
      npart_total   += hdr.npart;
      if(hdr.max_npart > max_npart)
	max_npart = hdr.max_npart;
      fprintf(stderr,"file = `%s' ngroup_total = %"PRId64"\n",argv[i],hdr.ngroups_total);
      fclose(fp);
    }  
  } else {
    ngroups_total = hdr.ngroups_total;
    npart_total = hdr.npart_total;
    max_npart = hdr.max_npart_total;
  }
  assert(ngroups_total > 0);
  
  hdr.num_files = 1 ;
  hdr.file_id = 0 ;
  hdr.format_part_data = 0 ;
  hdr.ngroups         = ngroups_total ;
  hdr.ngroups_total   = ngroups_total ;
  hdr.npart           = npart_total ;
  hdr.npart_total     = npart_total ;
  hdr.max_npart       = max_npart;
  hdr.max_npart_total = max_npart;
  hdr.bounds[0] = 0. ;
  hdr.bounds[1] = 0. ;
  hdr.bounds[2] = 0. ;
  hdr.bounds[3] = hdr.box_size ;
  hdr.bounds[4] = hdr.box_size ;
  hdr.bounds[5] = hdr.box_size ;
  
  bgc_write_header( stdout, hdr );

  uint32_t pad;
  pad = hdr.ngroups_total * bgc_sizeof_gdata( hdr.format_group_data ) ;
  bgc_write_raw_pad( stdout, pad) ;

  for(int i=1; i<argc; i++ ) {
    group_readwrite( argv[i] );
  }
  bgc_write_raw_pad( stdout, pad) ;


  fprintf(stderr,"bgc2_reduce> Done. Wrote %"PRId64" into groups zero_sat_count = %"PRId64"\n",hdr.ngroups_total,zero_sat_count);
  exit(EXIT_SUCCESS);
}
