
/* #include <stdlib.h> */
/* #include <stdio.h> */
/* #include <assert.h> */

//#include "binary_output.h"

#include "bgc_read_utils.h"
#include "utils.h"

static int BGC_VERBOSE = 0;

void bgc_read_header(FILE *fp, OUTPUT_HEADER * hdr)
{
    int pad;

    assert(fp != 0);

    fread(&pad,sizeof(int),1,fp);
    assert(pad == 1024);
    fread(hdr,sizeof(OUTPUT_HEADER),1,fp);
    fread(&pad,sizeof(int),1,fp);
    assert(pad == 1024);

    if(BGC_VERBOSE) {
        printf("READING HEADER INFORMATION:\n");
        printf("  total_files = %d\n", hdr->num_files);
        printf("  ngroups = %d\n", hdr->ngroups);
        printf("  starting at gid = %d\n", hdr->first_group_id);
        printf("  nparticles = %u\n", hdr->npart);
	printf("  nparticles_total = %u\n", hdr->npart_total);
	printf("  format     = %d\n", hdr->format);
        fflush(stdout);
    }
}

int * bgc_read_grouplist(FILE *fp, const OUTPUT_HEADER hdr)
{
    int pad;
    int *nParticlesPerGroup;
    /* int tmp; */
    //nParticlesPerGroup = calloc(hdr.ngroups+1, sizeof(int));
    nParticlesPerGroup = calloc(hdr.ngroups, sizeof(int));
    assert(nParticlesPerGroup != NULL);

    my_fread(&pad,sizeof(int),1,fp);
    assert(pad==sizeof(int)*hdr.ngroups);
    my_fread(nParticlesPerGroup,sizeof(int),hdr.ngroups,fp);
    /* for(i=0; i < hdr.ngroups; i++) { */
    /*     fread(&tmp,sizeof(int),1,fp); */
    /*     nParticlesPerGroup[i] = tmp; */
        /* if(BGC_VERBOSE) */
        /*     printf(" grp %4d: %5d\n", (i+hdr.first_group_id), tmp); */
    /* } */
    my_fread(&pad,sizeof(int),1,fp);
    assert(pad==sizeof(int)*hdr.ngroups);

    return nParticlesPerGroup;
}



/* Given a file pointer positioned after reading in the header, 
   skip the required number of bytes to reach the specified haloid.
   Need to read the header first -> know the number of groups+bytes
   to jump. haloid is the index in the file and not in the complete
   halo catalog (i.e., haloid does not include hdr.first_group_id
*/

void bgc_skip_to_haloid(FILE *fp, const OUTPUT_HEADER *hdr,int *nParticlesPerGroup,const unsigned int haloid)
{
  size_t bytes=0,size=0;
  int i;

  size = bgc_sizeof_pdata(hdr->format);
  
  for(i=0;i<haloid;i++) {
    bytes += sizeof(int);
    bytes += nParticlesPerGroup[i]*size;
    bytes += sizeof(int);
  }
  fseek(fp,bytes,SEEK_CUR);//done -> fp is now at the right location. The caller can verify

}


/* Read particle data for one group. One must cast the result appropriately */
void * bgc_read_particles(FILE *fp, const unsigned int npart, const int pdata_format)
{
    int pad;
    void *pd;
    size_t res;

    size_t size = bgc_sizeof_pdata(pdata_format);

    pd = malloc(npart * size);
    assert(pd != NULL);

    fread(&pad,sizeof(int),1,fp);
    assert(pad == npart * size);
    res = fread(pd,size,npart,fp);
    assert( res == npart );
    fread(&pad,sizeof(int),1,fp);
    assert(pad == npart * size);

    return (void*)pd;
}

/* Read particle data for one group. One must cast the result appropriately */
void bgc_read_part_into(FILE *fp, const unsigned int npart, const int pdata_format, void * pdata)
{
    int pad;
    size_t res;
    size_t size = bgc_sizeof_pdata(pdata_format);

    assert(pdata != NULL);
    
    fread(&pad,sizeof(int),1,fp);
    if( pad != npart * size ) {
      fprintf(stderr,"pad = %d (expected %zd)\n", pad, npart*size);
      fprintf(stderr,"npart = %d  particle_size = %zd\n", npart, size);
    }
    assert(pad == npart * size);
    res = fread(pdata,size,npart,fp);
    assert( res == npart );
    fread(&pad,sizeof(int),1,fp);
    assert(pad == npart * size);

    return;
}


//Not functional
/* void bgc_read_part_into_from_stream(char *buffer, const unsigned int npart, const int pdata_format, void * pdata) */
/* { */
/*     int pad; */
/*     /\* size_t res; *\/ */
/*     char int_as_char[5]; */
/*     size_t size = bgc_sizeof_pdata(pdata_format); */

/*     assert(pdata != NULL); */
    
    
/*     /\* fread(&pad,sizeof(int),1,fp); *\/ */
/*     memcpy(&pad,buffer,sizeof(int)); */
/*     /\* memcpy(int_as_char,buffer,sizeof(int)); *\/ */
/*     /\* int_as_char[4]='\0'; *\/ */
/*     /\* pad = strtol(int_as_char,NULL,10); *\/ */
/*     /\* fprintf(stderr,"buffer=%p\n",buffer); *\/ */
/*     /\* for(int i=0;i<5;i++) { *\/ */
/*     /\*   fprintf(stderr,"i=%d buf=%c \n",i,buffer[i]); *\/ */
/*     /\* } *\/ */

/*     /\* pad = *((int *) (buffer)); *\/ */
/*     if( pad != npart * size ) { */
/*       fprintf(stderr,"pad = %d (expected %zd)\n", pad, npart*size); */
/*       fprintf(stderr,"npart = %d  particle_size = %zd\n", npart, size); */
/*     } */
/*     assert(pad == npart * size); */
/*     buffer += sizeof(int); */
/*     //res = fread(pdata,size,npart,fp); */
/*     //assert( res == npart ); */
/*     memcpy(pdata,buffer,size*npart); */

/*     buffer += size*npart; */
/*     //fread(&pad,sizeof(int),1,fp); */
/*     memcpy(int_as_char,buffer,sizeof(int)); */
/*     int_as_char[4]='\0'; */
/*     pad = strtol(int_as_char,NULL,10); */

/*     if( pad != npart * size ) { */
/*       fprintf(stderr,"pad = %d (expected %zd)\n", pad, npart*size); */
/*       fprintf(stderr,"npart = %d  particle_size = %zd\n", npart, size); */
/*     } */

/*     assert(pad == npart * size); */
/*     buffer += sizeof(int); */
/*     return; */
/* } */


void print_pdata_format(FILE * fp, const int pdata_format)
{
    switch( pdata_format )
    {
        case PDATA_FORMAT_ID :
            fprintf(fp, "ID");
            break;
        case PDATA_FORMAT_IDBE :
            fprintf(fp,"IDBE");
            break;
        case PDATA_FORMAT_POS :
            fprintf(fp,"POS");
            break;
        case PDATA_FORMAT_POSBE :
            fprintf(fp,"POSBE");
            break;
        case PDATA_FORMAT_PV :
            fprintf(fp,"PV");
            break;
        case PDATA_FORMAT_PVBE :
            fprintf(fp,"PVBE");
            break;
        case PDATA_FORMAT_PVM :
            fprintf(fp,"PVM");
            break;
        case PDATA_FORMAT_PVMBE :
            fprintf(fp,"PVMBE");
            break;
        case PDATA_FORMAT_GPVM :
            fprintf(fp,"KITCHEN SINK (GPVM)");
            break;
        default :
            fprintf(stderr, "ERROR: unknown particle data format!  (format = %d)\n", pdata_format);
    }
}
