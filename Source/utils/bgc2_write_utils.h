#pragma once

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <stdint.h>
#include <inttypes.h>
#include <err.h>

#include "bgc2.h"
#include "bgc2_compat.h"

size_t ftwrite( void const *ptr, const size_t size, const size_t nitems, FILE * stream );
void bgc_write_header( FILE * fp, const OUTPUT_HEADER hdr );
void bgc_write_raw_pad( FILE * fp, const uint32_t pad );
void bgc_write_grouplist_raw( FILE * fp, const int64_t ngroups, const int64_t gdata_format, void *gd );
void bgc_write_grouplist( FILE * fp, const int64_t ngroups, const int64_t gdata_format, void *gd );
void bgc_write_pdata( FILE * fp, const int64_t npart, const int64_t pdata_format, void const *pdata );

