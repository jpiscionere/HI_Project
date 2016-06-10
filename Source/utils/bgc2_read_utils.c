#include "bgc2_read_utils.h"

static int BGC_VERBOSE = 0;     /* make this 1 to increase verbosity for debugging, etc. */

size_t
single_ftread( void *ptr, size_t size, size_t nitems, FILE * stream, int seek )
{
    uint32_t nbyte1, nbyte2;
    size_t nitem1;
    size_t res;

    char *msg_init = "ftread error";

    /* unformatted data FORMAT, typically 4-byte boundaries */
    res = fread( &nbyte1, sizeof( uint32_t ), 1, stream );
    if( res != 1 ) {
        fprintf( stderr, "%s: file empty? \n", msg_init );
        perror( msg_init );
        return ( res );
    }
    if( ( size_t ) nbyte1 > size * nitems ) {
        fprintf( stderr,
                 "%s: expected byte size does not match item*size nbyte = %" PRIu32 ", size = %"
                 PRIu64 ", nitems = %" PRIu64 "\n", msg_init, nbyte1, ( uint64_t ) size,
                 ( uint64_t ) nitems );
        return ( 0 );
    }

    assert( ( size > 0 ) && !( nbyte1 % size ) );

    if( seek ) {
        res = fseek( stream, size *  nbyte1, SEEK_CUR );
        nitem1 = ( res != 0 ) ? 0 : ( nbyte1 / size );
    } else {
        nitem1 = fread( ptr, size, nbyte1 / size, stream );
    }

    if( nitem1 != ( nbyte1 / size ) ) {
        fprintf( stderr,
                 "%s: %" PRIu64 " items expected, %" PRIu64 " items read. \n",
                 msg_init, ( uint64_t ) ( nbyte1 / size ), ( uint64_t ) nitem1 );
    }
    res = fread( &nbyte2, sizeof( int ), 1, stream );
    if( res != 1 ) {
        fprintf( stderr, "%s: file too short? \n", msg_init );
        perror( msg_init );
        return ( res );
    }
    if( nbyte1 != nbyte2 ) {
        fprintf( stderr,
                 "%s: byte paddings do not match, nbyte1 = %" PRIu32 ", nbyte2 = %" PRIu32 " \n",
                 msg_init, nbyte1, nbyte2 );
        return ( 0 );
    }

    return ( size_t ) nitem1;
}

size_t
_ftread( void *ptr, size_t size, size_t nitems, FILE * stream, int seek )
{
    char *msg_init = "ftread error";

    int64_t nread = 0, res;

    while( nread < nitems ) {
      res = single_ftread( ptr + ( size * nread ), size, nitems - nread, stream, seek );
      if( res <= 0 ) {
	fprintf( stderr,
		 "%s: expected byte size does not match item*size nbyte = %" PRIu64 ", size = %"
		 PRIu64 ", nitems = %" PRIu64 "\n", msg_init, ( uint64_t ) ( nread * size ),
		 ( uint64_t ) size, ( uint64_t ) nitems );
	return ( 0 );
      }
      nread += res;
    }
    return ( ( size_t ) nread );
}

/* Easy unformatted FORTRAN reading, with some sanity checking.
 * Usage like C fread. */
size_t
ftread( void *ptr, size_t size, size_t nitems, FILE * stream )
{
    return _ftread( ptr, size, nitems, stream, 0 );
}

void bgc_read_header( FILE * fp, OUTPUT_HEADER * hdr )
{
    int64_t size;

    assert( fp != 0 );
    size = ( int64_t ) ftread( hdr, sizeof( OUTPUT_HEADER ), 1, fp );
    assert( size == 1 );
    assert( hdr->magic == BGC_MAGIC );
    assert( hdr->version == 2 );

    if( BGC_VERBOSE ) {
        printf( "READING HEADER INFORMATION:\n" );
        printf( "  total_files = %" PRId64 "\n", hdr->num_files );
        printf( "  ngroups = %" PRId64 "\n", hdr->ngroups );
        printf( "  nparticles = %" PRId64 "\n", hdr->npart );
        printf( "  group data size = %zu\n", bgc_sizeof_gdata( hdr->format_group_data ) );
        fflush( stdout );
    }
}

void *
bgc_read_grouplist( FILE * fp, const OUTPUT_HEADER hdr )
{
    int64_t i;

    size_t size = bgc_sizeof_gdata( hdr.format_group_data );

    void *gd;

    GROUP_DATA_ID *halo;

    gd = calloc( hdr.ngroups + 1, size );
    assert( gd != NULL );

    ftread( gd, size, hdr.ngroups, fp );

    if( BGC_VERBOSE && ( size > 0 ) )
        for( i = 0; i < hdr.ngroups; i++ ) {
            halo = gd + i * size;
            printf( " grp %4" PRId64 ": %5" PRId64 "\n", halo->id, halo->npart );
        }

    return gd;
}

/* Read particle data for one group. One must cast the result appropriately */
void *
bgc_read_particles( FILE * fp, const uint64_t npart, const int64_t pdata_format )
{
    void *pd;

    size_t size = bgc_sizeof_pdata( pdata_format );

    pd = calloc( npart, size );
    assert( pd != NULL );

    ftread( pd, size, npart, fp );

    return ( void * )pd;
}

/* Read particle data for one group into *pdata.  Assumes memory is properly allocated to
 * contain FULL data.  One must cast the result appropriately to use it */
void
bgc_read_part_into( FILE * fp, const uint64_t npart, const int64_t pdata_format, void *pdata )
{
    size_t size = bgc_sizeof_pdata( pdata_format );

    assert( pdata != NULL );
    ftread( pdata, size, npart, fp );

    return;
}

/* skip over group without reading it */
int
bgc_skip_particles( FILE * fp, const uint64_t npart, const int64_t pdata_format )
{
    size_t size = bgc_sizeof_pdata( pdata_format ), res;

    res = _ftread( NULL, size, npart, fp, 1 );
    return ( ( res == npart ) ? 0 : -1 );
}

char * gdata_format_name( const int64_t gdata_format )
{
    switch ( gdata_format ) {
        case GDATA_FORMAT_ID:
            return "ID";
        case GDATA_FORMAT_RM:
            return "RM";
        case GDATA_FORMAT_RMPV:
            return "RMPV";
        case GDATA_FORMAT_RMPVMAX:
            return "RMPVMAX";
    }
    fprintf( stderr, "ERROR: unknown particle data format!  (format = %" PRId64 ")\n",
             gdata_format );
    return "";
}

void
print_format_group_data( FILE * fp, const int64_t gdata_format )
{
    fprintf( fp, "%s", gdata_format_name( gdata_format ) );
}

void
print_format_part_data( FILE * fp, const int64_t pdata_format )
{
    switch ( pdata_format ) {
        case PDATA_FORMAT_ID:
            fprintf( fp, "ID" );
            break;
        case PDATA_FORMAT_IDBE:
            fprintf( fp, "IDBE" );
            break;
        case PDATA_FORMAT_POS:
            fprintf( fp, "POS" );
            break;
        case PDATA_FORMAT_POSBE:
            fprintf( fp, "POSBE" );
            break;
        case PDATA_FORMAT_PV:
            fprintf( fp, "PV" );
            break;
        case PDATA_FORMAT_PVBE:
            fprintf( fp, "PVBE" );
            break;
        case PDATA_FORMAT_PVM:
            fprintf( fp, "PVM" );
            break;
        case PDATA_FORMAT_PVMBE:
            fprintf( fp, "PVMBE" );
            break;
        case PDATA_FORMAT_GPVM:
            fprintf( fp, "KITCHEN SINK (GPVM)" );
            break;
        default:
            fprintf( stderr, "ERROR: unknown particle data format!  (format = %" PRId64 ")\n",
                     pdata_format );
    }
}
