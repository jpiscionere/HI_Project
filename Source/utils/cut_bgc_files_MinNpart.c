#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "utils.h"
#include "bgc_read_utils.h"
#include "bgc_write_utils.h"


#define MAXLEN   (1000)

const char base_data_dir[]="/data0/LasDamas";
const char base_output_dir[]="/data0/LasDamas";

int main(int argc, char **argv)
{

  int lum_sample,seed,min_numpart;
  char fname[MAXLEN],outfname[MAXLEN],bgc_file_first_section[MAXLEN];
  FILE *fp=NULL,*outfp=NULL;
  OUTPUT_HEADER hdr,outheader ;
  int *nParticlesPerGroup = NULL;
  int **nParticlesPerGroup_local = NULL;
  PARTICLE_DATA_PV *pdata,pd ;
  int nfiles=0,ngroups;
  unsigned int *npart_local=NULL,npart_total;
  int *ngroups_local=NULL,ngroups_total;//number of particles, groups in the file
  int first_group_id=0;
  
  if(argc <= 2) { 
    fprintf(stderr,"ERROR: usage `%s'  <MinNpart ([1,99999]) > <list of 0000.bgc files>\n",argv[0]);
    for(int i=1;i<argc;i++)
      fprintf(stderr,"\t\t argv[%d] = %s \n",i,argv[i]);

    exit(EXIT_FAILURE);
  } else {
    min_numpart=atoi(argv[1]);
  }
  assert(min_numpart > 0);
  int max_ngroups=0;
  size_t size;
  for(int ifile=2;ifile<argc;ifile++) {
    my_snprintf(fname,MAXLEN,"%s",argv[ifile]);
    //my_snprintf(bgc_file_first_section,MAXLEN,"%*s",strlen(fname)-10,fname);
    size=strlen(fname)-9;
    strncpy(bgc_file_first_section,fname,size);
    bgc_file_first_section[size] = '\0';
    fprintf(stderr,"fname=%s strlen(fname)=%zu\n",fname,strlen(fname));
    fprintf(stderr,"bgc__=%s strlen(bgc_file_first_section) = %zu\n",bgc_file_first_section,strlen(bgc_file_first_section));
    fp=my_fopen(fname,"r");
    bgc_read_header(fp,&hdr);
    fclose(fp);
    exit(EXIT_SUCCESS);

    nfiles=hdr.num_files;
    pdata = my_malloc(bgc_sizeof_pdata(hdr.format), hdr.max_npart_total );

    fprintf(stderr,"mallocing for %d elements (nfiles) of ",nfiles);
    fprintf(stderr," %zu bytes\n",sizeof(*ngroups_local));
    ngroups_local = my_calloc(sizeof(*ngroups_local),nfiles);
    npart_local  = my_calloc(sizeof(*npart_local),nfiles);

    for(int i=0;i<nfiles;i++) {
      my_snprintf(fname,MAXLEN,"%s.%04d.bgc",bgc_file_first_section,i);
      fp=my_fopen(fname,"r");
      bgc_read_header(fp,&hdr);
      fclose(fp);
      ngroups=hdr.ngroups;
      if (ngroups > max_ngroups)
	max_ngroups=ngroups;
    }
    nParticlesPerGroup_local = (int **) matrix_malloc(sizeof(int),nfiles,max_ngroups);

    npart_total=0;
    ngroups_total=0;
    for(int i=0;i<nfiles;i++) {
      my_snprintf(fname,MAXLEN,"%s.%04d.bgc",bgc_file_first_section,i);
      fprintf(stderr,"Reading file `%s' ",fname);
      fp=my_fopen(fname,"r");
      bgc_read_header(fp,&hdr);
      ngroups=hdr.ngroups;
      fprintf(stderr,"ngroups = %d \n",ngroups);
      nParticlesPerGroup = bgc_read_grouplist(fp,hdr) ;      
      fclose(fp);
      ngroups_local[i]=0;
      for(int j=0;j<ngroups;j++) {
	if(nParticlesPerGroup[j] >=min_numpart) {
	  npart_local[i] += nParticlesPerGroup[j];
	  nParticlesPerGroup_local[i][ngroups_local[i]] = nParticlesPerGroup[j];
	  ngroups_local[i]++;
	}
      }
      ngroups_total += ngroups_local[i];
      npart_total += npart_local[i];
      free(nParticlesPerGroup);
    }

    fprintf(stderr,"Original (npart_total,ngroups_total) = (%u,%d) now = (%u,%d) \n",
	    hdr.npart_total,hdr.ngroups_total,
	    npart_total,ngroups_total);
    
    //Now output the data
    first_group_id=1;
    for(int i=0;i<nfiles;i++) {
      fprintf(stderr,"Writing file # %3d out %3d files. Ngroups = %8d Npart_local = %9d ...",i+1,nfiles,ngroups_local[i],npart_local[i]);
      my_snprintf(fname,MAXLEN,"%s.%04d.bgc",bgc_file_first_section,i);
      fp=my_fopen(fname,"r");
      bgc_read_header(fp,&hdr);
      nParticlesPerGroup = bgc_read_grouplist(fp,hdr) ;

      //now prepare the output buffer
      memcpy(&outheader,&hdr,sizeof(OUTPUT_HEADER));

      //update fields in the header that may have changed
      outheader.first_group_id=first_group_id;
      outheader.ngroups = ngroups_local[i];
      outheader.ngroups_total = ngroups_total;

      outheader.min_group_part = min_numpart;
      outheader.npart = npart_local[i];
      outheader.npart_total = npart_total;
      //done with updating the header
      
      my_snprintf(outfname,MAXLEN,"%s_MinNpart_%05d.%04d.bgc",bgc_file_first_section,min_numpart,i);
      outfp=my_fopen(outfname,"w");
      bgc_write_header(outfp,outheader);

      //write grouplist
      bgc_write_grouplist(outfp,ngroups_local[i],&(nParticlesPerGroup_local[i][0]));

      
      for(int j=0;j<hdr.ngroups;j++) {
	if(nParticlesPerGroup[j] >= min_numpart) {
	  //read the data 
	  bgc_read_part_into(fp, nParticlesPerGroup[j], hdr.format, pdata) ;

	  //output the particle data
	  bgc_write_pdata(outfp,nParticlesPerGroup[j],hdr.format,pdata);
	} else {
	  //skip over the bytes for this group
	  size_t bytes=0;
	  bytes += sizeof(int);//padding at beginning 
	  bytes += nParticlesPerGroup[j]*bgc_sizeof_pdata(hdr.format);
	  bytes += sizeof(int);//padding at end
	  my_fseek(fp,bytes,SEEK_CUR);
	}
      }
      fclose(fp);
      fclose(outfp);
      fprintf(stderr,"..done\n");
      first_group_id += ngroups_local[i];
      free(nParticlesPerGroup);
    }
    
    matrix_free((void **) nParticlesPerGroup_local,nfiles);
    free(ngroups_local);
    free(npart_local);
    free(pdata);
  }  

  return EXIT_SUCCESS;
}
