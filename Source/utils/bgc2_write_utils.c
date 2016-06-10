#include "bgc2_write_utils.h"

/* Easy unformatted FORTRAN writing, with some sanity checking.
 * Usage like C fwrite. */
size_t ftwrite( void const *ptr, const size_t size, const size_t nitems, FILE * stream )
{
    uint32_t nbytes;

    size_t nitem1;

    size_t res;

    char *msg_init = "ftwrite error";

    nbytes = ( uint32_t ) ( size * nitems );

    /* unformatted data FORMAT, typically 4-byte boundaries */
    res = fwrite( &nbytes, sizeof( uint32_t ), 1, stream );
    if( res != 1 ) {
        fprintf( stderr, "%s: file open? \n", msg_init );
        perror( msg_init );
        return ( res );
    }
    nitem1 = fwrite( ptr, size, nitems, stream );
    if( nitem1 != nitems ) {
        fprintf( stderr, "%s: %zu items requested, %zu items written. \n", msg_init,
                 nitems, nitem1 );
    }
    res = fwrite( &nbytes, sizeof( uint32_t ), 1, stream );
    if( res != 1 ) {
        fprintf( stderr, "%s: write error on second byte label\n", msg_init );
        perror( msg_init );
        return ( res );
    }

    return ( size_t ) nitems;
}


void bgc_write_header( FILE * fp, const OUTPUT_HEADER hdr )
{
    size_t res;

    /* sanity checks */
    assert( sizeof( OUTPUT_HEADER ) == OUTPUT_HEADER_SIZE );
    assert( fp != NULL );

    /* force BGC header to be at beginning of file! */
    rewind( fp );

    res = ftwrite( &hdr, sizeof( OUTPUT_HEADER ), 1, fp );

    assert( res == 1 );
}

void bgc_write_raw_pad( FILE * fp, const uint32_t pad )
{
    size_t res;

    res = fwrite( &pad, sizeof( uint32_t ), 1, fp );

    assert( res == 1 );
}

void bgc_write_grouplist_raw( FILE * fp, const int64_t ngroups, const int64_t gdata_format, void *gd )
{
    size_t res, size = bgc_sizeof_gdata( gdata_format );

    assert( gd != NULL );

    res = fwrite( gd, size, ngroups, fp );

    assert( res == ngroups );
}

void bgc_write_grouplist( FILE * fp, const int64_t ngroups, const int64_t gdata_format, void *gd )
{
    size_t res, size = bgc_sizeof_gdata( gdata_format );

    assert( gd != NULL );

    res = ftwrite( gd, size, ngroups, fp );

    assert( res == ngroups );
}

void bgc_write_pdata( FILE * fp, const int64_t npart, const int64_t pdata_format,
		      void const *pdata )
{
    size_t size, res;

    assert( pdata != NULL );

    size = bgc_sizeof_pdata( pdata_format );

    res = ftwrite( pdata, size, npart, fp );

    assert( res == npart );

}
