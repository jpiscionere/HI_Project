#pragma once

#include "ftread.h"
#include "ftwrite.h"

#include <stdint.h>
#include <inttypes.h>

enum DataType {dt_int32_t, dt_int64_t, dt_float, dt_double};

//Contains the data-types for the various columns 
/* extern int ColumnTypes[]; */
/* extern char *ColumnNames[]; */
/* extern int Ncolumns; */

void print_ascii_column(const enum DataType enumcolumntype, const void *extra_info, const int index, FILE *fp);
size_t sizeof_column_data(const enum DataType enumcolumntype);
