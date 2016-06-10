#ifndef _LOADSNAPSHOT_H_
#define _LOADSNAPSHOT_H_

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include "gadget_headers.h"

#ifndef id64
#define id64  unsigned long
#endif

#ifndef MAXLEN
#define MAXLEN    1000
#endif

#ifndef DOUBLE_EPS
#define   DOUBLE_EPS    (1e-10)
#endif





struct io_header get_gadget_header(const char *fname);
int get_gadget_nfiles(const char *fname);
struct particle_data * loadsnapshot(const char *fname,struct io_header *header);
void reordering(struct particle_data *P,id64 *Id,int64_t N);
int64_t get_Numpart(struct io_header *header);
void loadsnapshot_arrays(const char *fname,struct io_header *header,float **x,float **y,float **z,float **vx,float **vy,float **vz,id64 **part_ids);


#endif
