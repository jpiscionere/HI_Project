#ifndef __FTWRITE_H
#define __FTWRITE_H
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

/* int ftwrite(void * ptr, size_t size,size_t nitems, FILE * stream); */
size_t ftwrite(void const *ptr, const size_t size, const size_t nitems, FILE * stream );

#endif
