#ifndef __BGC_WRITE_UTILS_H__
#define __BGC_WRITE_UTILS_H__

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "bgc.h"
#include "ftwrite.h"

void bgc_write_header( FILE * fp, const OUTPUT_HEADER hdr );
void bgc_write_grouplist( FILE * fp, const int ngroups, int const *nParticlesPerGroup );
void bgc_write_pdata( FILE * fp, const unsigned int npart, const int pdata_format,void const *pdata );



#endif
