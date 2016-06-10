#ifndef __BGC_READ_UTILS
#define __BGC_READ_UTILS

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

//include the definitions of the particle structures etc..
#include "bgc.h"

void bgc_read_header(FILE *fp, OUTPUT_HEADER * hdr);
int * bgc_read_grouplist(FILE *fp, const OUTPUT_HEADER hdr);
void bgc_skip_to_haloid(FILE *fp, const OUTPUT_HEADER *hdr,int *nParticlesPerGroup,const unsigned int haloid);
void * bgc_read_particles(FILE *fp, const unsigned int npart, const int pdata_format);
inline void bgc_read_part_into(FILE *fp, const unsigned int npart, const int pdata_format, void * pdata);
/* inline void bgc_read_part_into_from_stream(char *buffer, const unsigned int npart, const int pdata_format, void * pdata); */
void print_pdata_format(FILE * fp, const int pdata_format);


#endif
