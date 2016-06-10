#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>

FILE
*check_fopen( const char *path, const char *mode )
{
    FILE *fp;

    fp = fopen(path, mode);

    if ( NULL == fp ) {
        char * msg;
        size_t len;

        len = strlen( path ) + strlen( mode ) + 40;

        msg = malloc( len * sizeof(char) );

        if( NULL == msg ) {
            fprintf( stderr, "ERROR: cannot allocate a string of %zd elements\n", len);
            exit( EXIT_FAILURE );
        }
        snprintf( msg , len, "could not open '%s' (mode %s)", path, mode );
        if( 0 == errno ) {
            fprintf( stderr, "ERROR: %s\n", msg );
        } else {
            perror( msg );
        }
        exit( EXIT_FAILURE );
    }

    return fp;
}

