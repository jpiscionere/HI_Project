#include "bgc_write_utils.h"


void bgc_write_header( FILE * fp, const OUTPUT_HEADER hdr )
{
  size_t res;

  /* sanity checks */
  assert( sizeof( OUTPUT_HEADER ) == OUTPUT_HEADER_SIZE );
  assert( fp != NULL );

  /* force BGC header to be at beginning of file! */
  rewind( fp );

  res = ftwrite( &hdr, sizeof( OUTPUT_HEADER ), 1, fp );

  assert( res == OUTPUT_HEADER_SIZE );
}

void bgc_write_grouplist( FILE * fp, const int ngroups, int const *nParticlesPerGroup )
{
  size_t res;

  assert( nParticlesPerGroup != NULL );

  res = ftwrite( nParticlesPerGroup, sizeof( int ), ngroups, fp );
  assert( res == sizeof( int ) * ngroups );
}

void bgc_write_pdata( FILE * fp, const unsigned int npart, const int pdata_format,
		      void const *pdata )
{
  size_t size, res;

  assert( pdata != NULL );

  size = bgc_sizeof_pdata( pdata_format );

  res = ftwrite( pdata, size, npart, fp );

  assert( res == size * npart );

}
