#include "loadsnapshot.h"
#include "utils.h"
#include "progressbar.h"
#include "sglib.h"
#include <omp.h>


// Assumes Gadget icformat=1.
struct io_header get_gadget_header(const char *fname)
{
  FILE *fp=NULL;
  char buf[MAXLEN], buf1[MAXLEN];
  int dummy;
  struct io_header header;
  my_snprintf(buf,MAXLEN, "%s.%d", fname, 0);
  my_snprintf(buf1,MAXLEN, "%s", fname);
  fp = fopen(buf,"r");
  if(fp == NULL) {
    fp = fopen(buf1,"r");
    if(fp == NULL) {
      fprintf(stderr,"ERROR: Could not find snapshot file.\n neither as `%s' nor as `%s'\n",buf,buf1);
      fprintf(stderr,"exiting..\n");
      exit(EXIT_FAILURE);
    }
  }

  //// Don't really care which file actually succeeded (as long as one, buf or buf1, is present)
  fread(&dummy, sizeof(dummy), 1, fp);
  fread(&header, sizeof(header), 1, fp);
  fread(&dummy, sizeof(dummy), 1, fp);
  fclose(fp);
  return header;
}

int get_gadget_nfiles(const char *fname)
{
  FILE *fp=NULL;
  char buf[MAXLEN], buf1[MAXLEN];
  int dummy;
  struct io_header header;
  my_snprintf(buf,MAXLEN, "%s.%d", fname, 0);
  my_snprintf(buf1,MAXLEN, "%s", fname);

  fp = fopen(buf,"r");
  if(fp == NULL) {
    fp = fopen(buf1,"r");
    if(fp == NULL) {
      fprintf(stderr,"ERROR: Could not find snapshot file.\n neither as `%s' nor as `%s'\n",buf,buf1);
      fprintf(stderr,"exiting..\n");
      exit(EXIT_FAILURE);
    }
  }

  //// Don't really care which file actually succeeded (as long as one, buf or buf1, is present)
  fread(&dummy, sizeof(dummy), 1, fp);
  fread(&header, sizeof(header), 1, fp);
  fread(&dummy, sizeof(dummy), 1, fp);
  fclose(fp);
  return header.num_files;
}


int64_t get_Numpart(struct io_header *header)
{
  int64_t N=0;
  if(header->num_files <= 1)
    for(int i = 0; i < 6; i++)
      header->npartTotal[i] = header->npart[i];

  for(int i = 0; i < 6; i++) {
    N += header->npartTotal[i];
    N += (((int64_t) header->npartTotalHighWord[i]) << 32);
  }

  return N;
}



/* this routine loads particle data from Gadget's default
 * binary file format. (A snapshot may be distributed
 * into multiple files.
 */

struct particle_data * loadsnapshot(const char *fname,struct io_header *header)
{
  FILE *fd=NULL;
  char buf[MAXLEN];
  int dummy,i;
  int64_t n,k,pc,pc_new;
  int64_t    ntot_withmasses;
  id64 *Id=NULL;
  size_t nmemb=0;
  const size_t OneElement=1;
  int files = get_gadget_nfiles(fname);
  int64_t NumPart = get_Numpart(header);
  struct particle_data *P=NULL;

#define SKIP my_fread(&dummy, sizeof(dummy),OneElement , fd);

  for(i=0, pc=0; i<files; i++, pc=pc_new) {
    if(files>1) {
      my_snprintf(buf,MAXLEN,"%s.%d",fname,i);
    } else {
      my_snprintf(buf,MAXLEN,"%s",fname);
    }

    fd = my_fopen(buf,"r");
    fprintf(stderr,"reading `%s' ...\n",buf); 
    
    nmemb = 1;
    my_fread(&dummy, sizeof(dummy), nmemb, fd);
    my_fread(header, sizeof(struct io_header), nmemb, fd);
    my_fread(&dummy, sizeof(dummy), nmemb, fd);
    
    for(k=0, ntot_withmasses=0; k<5; k++) {
      if(fabs(header->mass[k]) < DOUBLE_EPS) 
	ntot_withmasses+= header->npart[k];
    }

    if(i==0) {
      P  = my_malloc(sizeof(struct particle_data),NumPart);
      Id = my_malloc(sizeof(id64),NumPart);
    }

    SKIP;
    for(k=0,pc_new=pc;k<6;k++) {
      for(n=0;n<header->npart[k];n++) {
	nmemb = 3;
	my_fread(&P[pc_new].Pos[0], sizeof(float), nmemb, fd);
	pc_new++;
      }
    }
    SKIP;

    SKIP;
    for(k=0,pc_new=pc;k<6;k++) {
      for(n=0;n<header->npart[k];n++) {
	nmemb = 3;
	my_fread(&P[pc_new].Vel[0], sizeof(float), nmemb, fd);
	pc_new++;
      }
    }
    SKIP;
    
    
    SKIP;
    for(k=0,pc_new=pc;k<6;k++) {
      for(n=0;n<header->npart[k];n++) {
	my_fread(&Id[pc_new], sizeof(id64), OneElement, fd);
	pc_new++;
      }
    }
    SKIP;


    if(ntot_withmasses>0)
      SKIP;

    for(k=0, pc_new=pc; k<6; k++) {
      for(n=0;n<header->npart[k];n++) {
	P[pc_new].Type=k;
	
	if(fabs(header->mass[k]) < DOUBLE_EPS )
	  my_fread(&P[pc_new].Mass, sizeof(float), OneElement, fd);
	else
	  P[pc_new].Mass= header->mass[k];
	pc_new++;
      }
    }

    if(ntot_withmasses>0)
      SKIP;

    fclose(fd);
  }
  reordering(P,Id,NumPart);//so I can access by particle ids in assign_vxcm
  free(Id);
  return P;
}


void reordering(struct particle_data *P,id64 *Id,int64_t N)
{
  fprintf(stderr,"reordering....");
#define MULTIPLE_ARRAY_EXCHANGER(type,a,i,j) { SGLIB_ARRAY_ELEMENTS_EXCHANGER(struct particle_data,P,i,j); SGLIB_ARRAY_ELEMENTS_EXCHANGER(id64,Id,i,j); }
  SGLIB_ARRAY_QUICK_SORT(id64, Id, N,SGLIB_SAFE_NUMERIC_COMPARATOR,MULTIPLE_ARRAY_EXCHANGER);

  fprintf(stderr,"done.\n");
}




void loadsnapshot_arrays(const char *fname,struct io_header *header,float **x,float **y,float **z,float **vx,float **vy,float **vz,id64 **part_ids)
{
  FILE *fdpos=NULL,*fdvel=NULL,*fdids=NULL;
  char buf[MAXLEN];
  int dummy,i;
  int64_t n,k,pc,pc_new;
  int64_t ntot_withmasses;
  size_t nmemb=0;
  const size_t OneElement=1;
  int files = get_gadget_nfiles(fname);
  int64_t NumPart = get_Numpart(header);
  float tmp[3];
  int64_t npart_this_file=0;
  int interrupted;
  init_my_progressbar(NumPart,&interrupted);

  /* int64_t npart_per_file[files]; */
  int64_t cumulative_npart[files];
  memset(cumulative_npart,0,sizeof(int64_t)*files);

#ifdef SKIP
#undef SKIP  
#define SKIP my_fread(&dummy, sizeof(dummy),OneElement , fdpos);
#endif
  
  for(i=0; i<files; i++) {
    if(files>1) {
      my_snprintf(buf,MAXLEN,"%s.%d",fname,i);
    } else {
      my_snprintf(buf,MAXLEN,"%s",fname);
    }
    
    fdpos = my_fopen(buf,"r");
    /* fprintf(stderr,"\n reading header from `%s' ...\n",buf);  */

    nmemb = 1;
    my_fread(&dummy, sizeof(dummy), nmemb, fdpos);
    my_fread(header, sizeof(struct io_header), nmemb, fdpos);
    my_fread(&dummy, sizeof(dummy), nmemb, fdpos);

    npart_this_file=0;
    for(k=0,ntot_withmasses=0; k<6; k++) {
      if(fabs(header->mass[k]) < DOUBLE_EPS) 
	ntot_withmasses+= header->npart[k];
      
      npart_this_file += header->npart[k];
    }
    assert(ntot_withmasses == 0);

    /* npart_per_file[i] = npart_this_file; */
    cumulative_npart[i] = npart_this_file;
    fclose(fdpos);
  }

  //get the cumulative particle load
  for(i=1;i<files;i++) {
    cumulative_npart[i] += cumulative_npart[i-1];
  }

  *x  = my_malloc(sizeof(float),NumPart);
  *y  = my_malloc(sizeof(float),NumPart);
  *z  = my_malloc(sizeof(float),NumPart);
  *vx = my_malloc(sizeof(float),NumPart);
  *vy = my_malloc(sizeof(float),NumPart);
  *vz = my_malloc(sizeof(float),NumPart);
  *part_ids = my_malloc(sizeof(id64),NumPart);

  /* id64 xx;   */
  size_t bytes_to_skip=0;

/* #pragma omp parallel num_threads(files)  private(i,pc,pc_new,buf,dummy,fdpos,header,npart_this_file,tmp,k,nmemb,xx,fdvel,fdids) */
/*   { */
/*     int tid=omp_get_thread_num(); */

/*     if(tid < files) { */
  {
    for(i=0;i<files;i++) {
      /* i = tid; */
      if(files>1) {
	my_snprintf(buf,MAXLEN,"%s.%d",fname,i);
      } else {
	my_snprintf(buf,MAXLEN,"%s",fname);
      }
      if(files>1) {
	my_snprintf(buf,MAXLEN,"%s.%d",fname,i);
      } else {
	my_snprintf(buf,MAXLEN,"%s",fname);
      }
      /* fprintf(stderr,"tid=%3d buf = %s \n",tid,buf);       */
      if(i==0) {
	pc = 0;
      } else {
	pc = cumulative_npart[i-1];
      }
      
      fdpos = my_fopen(buf,"r");
      /* fprintf(stderr,"\nreading `%s' ... pc = %15"PRId64"\n",buf,pc); */
      /* interrupted=1; */

      //Reading begins
      nmemb = 1;
      my_fread(&dummy, sizeof(dummy), nmemb, fdpos);
      my_fread(header, sizeof(struct io_header), nmemb, fdpos);
      my_fread(&dummy, sizeof(dummy), nmemb, fdpos);
    
      npart_this_file=0;
      for(k=0; k<6; k++) {
	npart_this_file += header->npart[k];
      }
    
      fdvel = my_fopen(buf,"r");
      fdids = my_fopen(buf,"r");

      //now point the fdvel + fdids to the right places

      //Position the velocity file pointer
      bytes_to_skip = 4 + 256 +  4;//Takes care of the header
      bytes_to_skip += 4;//for the pad for the positions
      bytes_to_skip += npart_this_file*sizeof(float)*3 ;
      bytes_to_skip += 4;//for the pad for the positions
    
      /* fprintf(stderr,"npart_this_fil = %"PRId64" bytes_to_skip = %zu \n",npart_this_file,bytes_to_skip); */
    
      my_fseek(fdvel,(long)bytes_to_skip,SEEK_CUR);
      nmemb=1;
      my_fread(&dummy,sizeof(dummy),nmemb,fdvel);
      /* fprintf(stderr,"dummy = %d 3*sizeof(float)*npart_this_file = %zu\n",dummy,3*sizeof(float)*npart_this_file); */
      assert(dummy == 3*sizeof(float)*npart_this_file);


      //Now position the id file pointer
      bytes_to_skip +=4; //for the initial padding for the velocities
      bytes_to_skip += npart_this_file*sizeof(float)*3 ;
      bytes_to_skip += 4;//for the final padding for the velocities
      my_fseek(fdids,(long)bytes_to_skip,SEEK_CUR);

      nmemb=1;
      my_fread(&dummy,sizeof(dummy),nmemb,fdids);
      assert(dummy == sizeof(id64)*npart_this_file);
    
    
      SKIP;
      assert(dummy == 3*sizeof(float)*npart_this_file);
      pc_new = pc;
      for(k=0;k<6;k++) {
	for(n=0;n<header->npart[k];n++) {
	  my_progressbar(pc_new,&interrupted);
	  nmemb = 3;
	  my_fread(tmp, sizeof(float), nmemb, fdpos);
	  /* for(int ii=0;ii<3;ii++) { */
	  /*   assert(tmp[ii] >=0.0 && tmp[ii]  <= header->BoxSize); */
	  /* } */

	  (*x)[pc_new] = tmp[0];
	  (*y)[pc_new] = tmp[1];
	  (*z)[pc_new] = tmp[2];

	  my_fread(tmp, sizeof(float), nmemb, fdvel);
	  (*vx)[pc_new] = tmp[0];
	  (*vy)[pc_new] = tmp[1];
	  (*vz)[pc_new] = tmp[2];
		  pc_new++;
		}
      }
      my_fread((*part_ids), sizeof(id64), npart_this_file,fdids);
      
      fclose(fdpos);
      fclose(fdvel);
      fclose(fdids);
    }
    finish_myprogressbar(&interrupted);

  }
}




