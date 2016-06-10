#include "fastfood_extra.h"

/* int  ColumnTypes[] = {dt_double, dt_int32_t, dt_int64_t};//valid types are dt_double or dt_float, and dt_int32_t and dt_int64_t */
/* char *ColumnNames[] = {"halomass","central_flag","haloid"}; */
/* int  Ncolumns = sizeof(ColumnTypes)/sizeof(ColumnTypes[0]); */

void print_ascii_column(const enum DataType enumcolumntype, const void *extra_info, const int index, FILE *fp)
{
  assert(extra_info != NULL);
  switch ( enumcolumntype ) {
  case dt_double:
    {
      const double *tmp = extra_info;
      fprintf(fp, " %e ", tmp[index]);
      break;
    }
  case dt_int64_t:
    {
      const int64_t *tmp = extra_info;
      fprintf(fp," %"PRId64" ",tmp[index]);
      break;
    }
  case dt_float:
    {
      const float *tmp = extra_info;
      fprintf(fp," %e ",tmp[index]);
      break;
    }
  case dt_int32_t:
    {
      const int32_t *tmp = extra_info;
      fprintf(fp," %d ",tmp[index]);     
      break;
    }
  default:
    {
      fprintf(stderr,"WARNING: ColumnType %d is unknown \n",enumcolumntype);
      break;
    }
  }
}


size_t sizeof_column_data(const enum DataType enumcolumntype)
{

  switch ( enumcolumntype) {
  case dt_double:
    return sizeof(double);
  case dt_int64_t:
    return sizeof(int64_t);
  case dt_float:
    return sizeof(float);
  case dt_int32_t:
    return sizeof(int32_t);
  }

  fprintf( stderr, "ERROR: unknown column type !  (format = %d\n", enumcolumntype);
  return 0;
}


