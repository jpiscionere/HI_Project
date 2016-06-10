/*   	ftwrite writes data unformatted using fortran convention --
	i.e. an integer specifying the number of bytes in the record,
	the data record, and another integer specifying the number of 
	bytes in the record.  The call is identical to the standard
	i/o library routine fwrite.
*/

#include "ftwrite.h"

/* int ftwrite(void * ptr, size_t size,size_t nitems, FILE * stream) */
/* { */
/*   int nbytes; */
/*   size_t nitem1 ; */
/*   int errno ; */
  
/*   errno = 0 ; */
/*   nbytes = size*nitems ; */
/*   if ( fwrite(&nbytes,sizeof(int),1,stream) != 1 )  { */
/*     errno = -10 ; */
/*     fprintf(stderr,"write error, is the file open ? \n") ; */
/*   } */
/*   nitem1 = fwrite(ptr,size,nitems,stream) ; */
/*   if ( nitem1 != nitems ) { */
/*     errno = -20 ; */
/*     fprintf(stderr,"write error, %zu items requested, %zu items written. \n", */
/* 	    nitems,nitem1) ; */
/*   } */
/*   if ( fwrite(&nbytes,sizeof(int),1,stream) != 1 )  { */
/*     errno = -30 ; */
/*     fprintf(stderr,"write error on second byte label \n") ; */
/*   } */
  
/*   return(errno) ; */
/* } */


/* Easy unformatted FORTRAN reading, with some sanity checking.
 * Usage like C fread. but returns size of main object read in or "0" on an error */
size_t ftwrite(void const *ptr, const size_t size, const size_t nitems, FILE * stream )
{
  int nbytes;
  size_t nitem1;
  size_t res;
  char *msg_init = "ftwrite error";

  nbytes = ( int )size *nitems;

  assert( nbytes > 0 );

  /* unformatted data FORMAT, typically 4-byte boundaries */
  res = fwrite( &nbytes, sizeof( int ), 1, stream );
  if( res != 1 ) {
    fprintf( stderr, "%s: file open? \n", msg_init );
    perror( msg_init );
    return ( res );
  }
  nitem1 = fwrite( ptr, size, nitems, stream );
  if( nitem1 != nitems ) {
    fprintf( stderr, "%s: %zu items requested, %zu items written. \n", msg_init, nitems,
	     nitem1 );
  }
  res = fwrite( &nbytes, sizeof( int ), 1, stream );
  if( res != 1 ) {
    fprintf( stderr, "%s: write error on second byte label\n", msg_init );
    perror( msg_init );
    return ( res );
  }

  return ( size_t ) nbytes;
}
