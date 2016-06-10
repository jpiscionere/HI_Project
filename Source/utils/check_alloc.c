#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <stdarg.h>

int my_snprintf(char *buffer,int len,const char *format, ...)
{
  va_list args;
  int nwritten=0;

  va_start(args,format);
  nwritten=vsnprintf(buffer, len, format, args );
  va_end(args);
  if (nwritten > len || nwritten < 0)
    {
      fprintf(stderr,"ERROR: printing to string failed (wrote %d characters while only %d characters were allocated)\n",nwritten,len);
      fprintf(stderr,"Increase maxlen in `defs.h' ..exiting\n");
      exit(EXIT_FAILURE);
    }
  return nwritten;
}



static inline
void *
check_realloc( void * data, size_t count, size_t size )
{
    /* first check to see what function we want */
    if( NULL == data ) {
        data = calloc( count, size );
    } else {
        data = realloc( data, count * size );
    }

    /* verify that the allocation worked */
    if( NULL == data ) {
        fprintf(stderr,
                "ERROR: could not allocate memory for %zu elements of %zu size\n",
                count, size);
        exit(1);
    }

    return data;
}

static inline
void *
check_alloc( size_t count, size_t size )
{
    void * data;

    data = check_realloc(NULL, count, size);

    return data;
}

