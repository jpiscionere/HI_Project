#pragma once

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <assert.h>
#include <inttypes.h>

#include "bgc2.h"
#include "bgc2_compat.h"

size_t single_ftread( void *ptr, size_t size, size_t nitems, FILE * stream, int seek );
size_t _ftread( void *ptr, size_t size, size_t nitems, FILE * stream, int seek );
size_t ftread( void *ptr, size_t size, size_t nitems, FILE * stream );
void bgc_read_header( FILE * fp, OUTPUT_HEADER * hdr );
void * bgc_read_grouplist( FILE * fp, const OUTPUT_HEADER hdr );
void * bgc_read_particles( FILE * fp, const uint64_t npart, const int64_t pdata_format );

void bgc_read_part_into( FILE * fp, const uint64_t npart, const int64_t pdata_format, void *pdata );
int bgc_skip_particles( FILE * fp, const uint64_t npart, const int64_t pdata_format );
char * gdata_format_name( const int64_t gdata_format );
void print_format_group_data( FILE * fp, const int64_t gdata_format );
void print_format_part_data( FILE * fp, const int64_t pdata_format );














