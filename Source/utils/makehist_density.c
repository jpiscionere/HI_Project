/* Code to generate an uniform grid of number density 
out of a Gadget snapshot (multiple particle types, 
multiple files per snapshot). The only thing you 
need to set is ID_BYTES in main().

Compiling (there should not be any warnings):
gcc -m64 -O3 -std=c99 -g  -Wextra -Wall -Wformat=3  -g -m64 -std=c99  -Wpacked  -Wnested-externs -Wpointer-arith  -Wredundant-decls  -Wfloat-equal -Wcast-qual  -Wshadow  -Wcast-align -Wmissing-declarations -Wmissing-prototypes  -Wnested-externs -fno-strict-aliasing -Wpadded -Wstrict-prototypes -Wconversion makehist_density.c -lm -o makehist_density



Inputs:
1. Nbins    -- number of bins. Output grid will be Nbins x Nbins x Nbins [integer]
2. inpath   -- input directory that contains the snapshots.              [string]
3. snapshot -- the snapshot number to analyze                            [integer]
4. basename -- the basename of snapshots                                 [string]
5. outpath  -- output directory (You need to have write permissions)     [string]

Output:
The output is a F77 unformatted binary file with the following structure: 

a) int Nbins                (padded with a size_t, 4 or 8 bytes)
b) struct gadget io_header  (contains the snapshot header from either snapshot or snaphot.0, same padding )
c) Nbins x Nbins x Nbins of doubles -> the number densities normalized to the total mass in the box.   

Bugs:
I have checked with valgrind, there are no memory leaks, un-initialized values etc. So
this code is fairly safe to use. In the off-chance you do actually find a bug,
please send it to me: manodeep at gmail dot com.

License:
Free to use, modify and distribute. So, my guess is GPL. 

Version 1.0 : Manodeep Sinha, 2nd April, 2012.
*/


/* The following #defines are to guard against 32 bit 
compilation and subsequent large file handling. You
should really upgrade to a 64 bit machine!! 

The lines are taken from Stackoverflow:
http://stackoverflow.com/questions/1505582/determining-32-vs-64-bit-in-c

I have intentionally disregarded checking for Windows. 

-MS 2nd April, 2012.
*/

//check 64 bit
#if __GNUC__
#if __x86_64__ || __ppc64__
#define ENVIRONMENT64
#else
#define ENVIRONMENT32  
#endif
#endif

//if 32 bit, enable large file macros
#ifdef ENVIRONMENT32
#define _LARGEFILE_SOURCE
#define _LARGEFILE64_SOURCE
#define _FILE_OFFSET_BITS 64
#endif


//regular code commences from here
#include<stdio.h>
#include<stdlib.h>
#include<stdint.h>
#include<inttypes.h>
#include<math.h>
#include<string.h>
#include<limits.h>
#include<time.h>
#include<stdarg.h>

const size_t MAXLEN=1000;//number of characters in a string (pathnames, filenames etc)
short FLAG_INDIVIDUAL_PARTMASS=0;//flag if the snapshot has individual particle masses.

struct io_header
{
  int32_t npart[6];                    /*!< number of particles of each type in this file */
  double mass[6];                      /*!< mass of particles of each type. If 0, then the masses are explicitly
										 stored in the mass-block of the snapshot file, otherwise they are omitted */
  double time;                         /*!< time of snapshot file */
  double redshift;                     /*!< redshift of snapshot file */
  int32_t flag_sfr;                    /*!< flags whether the simulation was including star formation */
  int32_t flag_feedback;               /*!< flags whether feedback was included (obsolete) */
  uint32_t npartTotal[6];              /*!< total number of particles of each type in this snapshot. This can be
										 different from npart if one is dealing with a multi-file snapshot. */
  int32_t flag_cooling;                /*!< flags whether cooling was included  */
  int32_t num_files;                   /*!< number of files in multi-file snapshot */
  double BoxSize;                      /*!< box-size of simulation in case periodic boundaries were used */
  double Omega0;                       /*!< matter density in units of critical density */
  double OmegaLambda;                  /*!< cosmological constant parameter */
  double HubbleParam;                  /*!< Hubble parameter in units of 100 km/sec/Mpc */
  int32_t flag_stellarage;             /*!< flags whether the file contains formation times of star particles */
  int32_t flag_metals;                 /*!< flags whether the file contains metallicity values for gas and star particles */
  uint32_t npartTotalHighWord[6];      /*!< High word of the total number of particles of each type */
  int32_t  flag_entropy_instead_u;     /*!< flags that IC-file contains entropy instead of u */
  char fill[60];                       /*!< fills to 256 Bytes */
};  


//the function declarations -- to avoid compile time warnings. 
//utilities

// MS: 08/14/2012 - check_string_copy is now obsolete.
/* inline void check_string_copy(const int nwritten,const size_t allocated_length); */
int my_snprintf(char *buffer,size_t len,const char *format, ...);
void my_fwrite(void *ptr, size_t size, size_t nmemb, FILE *stream);//note that these i/o functions are void as opposed to the standard ones that return size_t
void my_fread(void *ptr, size_t size, size_t nmemb, FILE *stream);
void * my_malloc(size_t size,int64_t N);
FILE * my_fopen(const char *fname,const char *mode);
void print_time(time_t t0,time_t t1,const char *s);

//functions specific to gadget snapshots
struct io_header get_gadget_header(const char *fname);
int64_t get_Numpart(struct io_header *header);


//Actual utilities being called. 
void my_fwrite(void *ptr, size_t size, size_t nmemb, FILE *stream)
{
  size_t nwritten;
  nwritten = fwrite(ptr, size, nmemb, stream);
  if(nwritten != nmemb)
    {
      fprintf(stderr,"I/O error (fwrite) has occured.\n");
	  fprintf(stderr,"Instead of reading nmemb=%zu, I got nread = %zu ..exiting\n",nmemb,nwritten);
      exit(EXIT_FAILURE);
    }
}


void my_fread(void *ptr, size_t size, size_t nmemb, FILE *stream)
{
  size_t nread;
  nread = fread(ptr, size, nmemb, stream);
  if(nread != nmemb)
    {
      fprintf(stderr,"I/O error (fread) has occured.\n");
	  fprintf(stderr,"Instead of reading nmemb=%zu, I got nread = %zu ..exiting\n",nmemb,nread);
      exit(EXIT_FAILURE);
    }
}



void * my_malloc(size_t size,int64_t N)
{
  size_t bytes = size*N;
  void *x = malloc(bytes);
  if(x==NULL)
	{
	  fprintf(stderr,"ERROR: Failed to allocate memory for %"PRId64" elements of size %zu bytes\n exiting\n",N,size);
	  exit(EXIT_FAILURE);
	}

  return x;
}

FILE * my_fopen(const char *fname,const char *mode)
{
  FILE *fp=NULL;
  fp = fopen(fname,mode);
  if(fp==NULL)
	{
	  fprintf(stderr,"Could not open file `%s'\n",fname);
	  exit(EXIT_FAILURE);
	}
  
  return fp;
}

struct io_header get_gadget_header(const char *fname)
{
  FILE *fp=NULL;
  char buf[MAXLEN], buf1[MAXLEN];
  int dummy;
  size_t one=1;
  struct io_header header;
  my_snprintf(buf,MAXLEN, "%s.%d", fname, 0);
  my_snprintf(buf1,MAXLEN, "%s", fname);

  if(sizeof(struct io_header) != 256)
	{
	  fprintf(stderr,"ERROR: Gadget header is not %zu bytes and not *exactly* 256 bytes..exiting\n",sizeof(struct io_header));
	  exit(EXIT_FAILURE);
	}

  fp = fopen(buf,"r");
  if(fp == NULL)
	{
	  fp = fopen(buf1,"r");
	  if(fp == NULL)
		{
		  fprintf(stderr,"ERROR: Could not find snapshot file.\n neither as `%s' nor as `%s'\n",buf,buf1);
		  fprintf(stderr,"exiting..\n");
		  exit(EXIT_FAILURE);
		}
	}

  //// Don't really care which file actually succeeded (as long as one, buf or buf1, is present)
  my_fread(&dummy, sizeof(dummy), one, fp);
  my_fread(&header, sizeof(header), one, fp);
  my_fread(&dummy, sizeof(dummy), one, fp);
  fclose(fp);
  return header;

}

int64_t get_Numpart(struct io_header *header)
{
  int64_t NumPart=0;

  if(header->num_files <= 1)
	for(int i = 0; i < 6; i++)
	  header->npartTotal[i] = header->npart[i];
  
  for(int i = 0; i < 6; i++)
	{
	  NumPart += header->npartTotal[i];
	  NumPart += (((int64_t) header->npartTotalHighWord[i]) << 32);

	  //check if there are particles with individual particles masses
	  if(header->npartTotal[i] > 0 && !(header->mass[i] > 0.0))
		FLAG_INDIVIDUAL_PARTMASS = 1; //set the flag indicating that particles have individual masses
	}

  return NumPart;
}


void print_time(time_t t0,time_t t1,const char *s)
{
  double timediff = difftime(t1,t0);
  double ratios[] = {24*3600.0,  3600.0,  60.0,  1};
  char units[4][10]  = {"days", "hrs" , "mins", "secs"};
  int which = 0;

  double timeleft = timediff;
  double time_to_print;
  fprintf(stderr,"Time taken to execute '%s'  = ",s);

  if(timediff < ratios[2])
	fprintf(stderr,"%5d secs",(int)timediff);
  else
	while (which < 4)
	  {
		time_to_print = floor(timeleft/ratios[which]);
		if (time_to_print > 1)
		  {
			timeleft -= (time_to_print*ratios[which]);
			fprintf(stderr,"%5d %s",(int)time_to_print,units[which]);
		  }
		which++;
	  }

  fprintf(stderr,"\n");

}


// A real wrapper to snprintf that will exit() if the allocated buffer length 
// was not sufficient. Usage is the same as snprintf 
int my_snprintf(char *buffer,size_t len,const char *format, ...)
{
  va_list args;
  int nwritten=0;
 
  va_start(args,format);
  nwritten=vsnprintf(buffer, len, format, args );
  va_end(args);  
  if ((size_t) nwritten > len || nwritten < 0)
	{
	  fprintf(stderr,"ERROR: printing to string failed (wrote %d characters while only %zu characters were allocated)\n",nwritten,len);
	  fprintf(stderr,"Increase MAXLEN ..exiting\n");
	  exit(EXIT_FAILURE);
	}
  return nwritten;
}


int main(int argc, char **argv)
{
  char inpath[MAXLEN],outpath[MAXLEN],snapshot_base[MAXLEN],snapshot_name[MAXLEN],outfname[MAXLEN];
  int snapshot,Nbins;

  char progressbarstring[MAXLEN];
  int interrupted=0;
  int64_t PRINTSTEP,SMALLPRINTSTEP,index;
  int percent,end_of_string_index;;

  //set this to 8 if the particle ID's are 8 bytes. 
  const size_t ID_BYTES = 8;


  int nfiles,ifile;
  struct io_header header,header1;
  int64_t NumPart;
  FILE *fd=NULL,*fdmass=NULL;
  double BoxSize,binsize;
  int dummy;
  float pos[3],partmass;

  double ***histogram=NULL;
  int xi,yi,zi;
  double total_part_mass=0.0;
  size_t one=1;
  time_t t_codestart,t_codeend,t_sectionstart,t_sectionend;

  t_codestart=time(NULL);
  if (argc !=6)
    {
      fprintf(stderr,"\n\n wrong argument(s).  Specify: `%s' \n\n",argv[0]);
      fprintf(stderr,"<Nbins>       (Nbins for the 3d histogram)\n");
      fprintf(stderr,"<inpath>      (input path without trailing slash)\n");
      fprintf(stderr,"<basename>    (basename of snapshot files)\n");
      fprintf(stderr,"<snapshot>    (number of snapshot)\n");
      fprintf(stderr,"<outpath>     (output path)\n");
      fprintf(stderr,"\n exiting ..\n");
      exit(EXIT_FAILURE);
    }
  else
    {
      Nbins = atoi(argv[1]);

	  //copy the input directory (containing the snapshots)
      my_snprintf(inpath,MAXLEN,"%s",argv[2]);

	  //copy the basename of the snapshots
	  my_snprintf(snapshot_base,MAXLEN,"%s",argv[3]);

	  //get the snapshot number
      snapshot=atoi(argv[4]);

	  //copy the output directory to outpath
	  my_snprintf(outpath,MAXLEN,"%s",argv[5]);

	  fprintf(stderr,"Running `%s' with the parameters:\n\n",argv[0]);
	  fprintf(stderr,"\t\t <Nbins>    = %d \n",Nbins);
	  fprintf(stderr,"\t\t <inpath>   = %s \n",inpath);
	  fprintf(stderr,"\t\t <basename> = %s \n",snapshot_base);
      fprintf(stderr,"\t\t <snapshot> = %d \n",snapshot);
	  fprintf(stderr,"\t\t <outpath>  = %s \n",outpath);
	  fprintf(stderr,"\n");
    }

  my_snprintf(snapshot_name,MAXLEN,"%s/%s_%03d",inpath,snapshot_base,snapshot);

  header  =  get_gadget_header(snapshot_name);
  nfiles  =  header.num_files;
  NumPart = get_Numpart(&header);
  if(NumPart > INT_MAX && ID_BYTES < 8)
	{
	  //logically, if NumPart > INT_MAX, ID should be of long long (8 bytes). 
	  //However, ID bytes can be 8 bytes even though NumPart is < INT_MAX. So, an "AND" 
	  //is required.
	  fprintf(stderr,"ERROR: Total number of particles = %"PRId64" exceeds INT_MAX = %d but ID bytes are set to %zu ..exiting\n",
			  NumPart,INT_MAX, ID_BYTES);
	  exit(EXIT_FAILURE);
	}

  if(NumPart > INT_MAX && sizeof(size_t) < 8)
	{
	  fprintf(stderr,"ERROR: Total number of particles = %"PRId64" exceeds INT_MAX = %d but the code has *NOT* been compiled in 64bit mode (sizeof(size_t) = %zu bytes) ..exiting\n",
			  NumPart,INT_MAX, sizeof(size_t));
	  fprintf(stderr,"[Hint: If you are using gcc, add the flag -m64 to force 64bit compilation]\n");
	  exit(EXIT_FAILURE);
	}

  BoxSize = header.BoxSize;
  binsize = BoxSize/Nbins;

  //allocate memory for histogram
  fprintf(stderr,"Allocating memory for the histogram ...");
  histogram = my_malloc(sizeof(histogram),(int64_t) Nbins);
  for(int i=0;i<Nbins;i++)
    {
      histogram[i] = my_malloc(sizeof(histogram[i]),(int64_t) Nbins);
      for(int j=0;j<Nbins;j++)
		{
		  histogram[i][j] = my_malloc(sizeof(histogram[i][j]),(int64_t) Nbins);
		  memset(histogram[i][j],0,sizeof(histogram[i][j])*Nbins);
		}
    }

  fprintf(stderr,"..done\n");

  //progressbarstring stores the entire string -> if 
  //the output is interrupted, then the whole progressbar
  //can be reprinted. 
  my_snprintf(progressbarstring,MAXLEN,"%s","0%");

  end_of_string_index = 2;
  interrupted = 0;

  fprintf(stderr,"\n\nCalculating the histogram for Gadget snapshots ... \n");
  PRINTSTEP = (int64_t) floor(0.1*NumPart);
  SMALLPRINTSTEP = ceil(0.01*NumPart) > 1 ? ceil(0.01*NumPart):1;

  ifile =0;
  index=0;
  while(ifile < nfiles)
	{
	  t_sectionstart = time(NULL);
	  if (nfiles == 1) //Ideally no need to check this string copy
		{
		  my_snprintf(snapshot_name,MAXLEN,"%s/%s_%03d",inpath,snapshot_base,snapshot);
		}
	  else
		{
		  if(nfiles > 1)
			{
			  my_snprintf(snapshot_name,MAXLEN,"%s/%s_%03d.%d",inpath,snapshot_base,snapshot,ifile);
			}
		  else
			{
			  fprintf(stderr,"number of files for a snapshot *HAS* to be at least 1. I got nfiles = %d .exiting\n",nfiles);
			  exit(EXIT_FAILURE);
			}
		}

	  fd = my_fopen(snapshot_name,"r");
	  fprintf(stderr,"Reading file `%s' ...",snapshot_name);
	  interrupted=1;
	  my_fread(&dummy, sizeof(dummy), one, fd);
	  my_fread(&header1, sizeof(header1), one, fd);
	  my_fread(&dummy, sizeof(dummy), one, fd);

	  my_fread(&dummy, sizeof(dummy), one, fd);//set it ready to read in the first of the particle positions
	  
	  if(FLAG_INDIVIDUAL_PARTMASS == 1)
		{
		  fdmass = my_fopen(snapshot_name,"r");
		  my_fread(&dummy, sizeof(dummy), one, fdmass);
		  my_fread(&header1, sizeof(header1), one, fdmass);
		  my_fread(&dummy, sizeof(dummy), one, fdmass);

		  /*skip position data */
		  my_fread(&dummy, sizeof(dummy), one, fdmass);
		  for(int k=0; k<6; k++) 
			fseek(fdmass, (long) (header1.npart[k]*sizeof(float)*3), SEEK_CUR);//casting to 'long' to avoid the signed/unsigned compilation warning.
		  my_fread(&dummy, sizeof(dummy), one, fdmass);

		  /* skip velocity data */
		  my_fread(&dummy, sizeof(dummy), one, fdmass);
		  for(int k=0; k<6; k++) 
			fseek(fdmass, (long) (header1.npart[k]*sizeof(float)*3), SEEK_CUR);//casting to 'long' to avoid the signed/unsigned compilation warning.
		  my_fread(&dummy, sizeof(dummy), one, fdmass);

		  /* skip id data */
		  my_fread(&dummy, sizeof(dummy), one, fdmass);
		  //different scope to check for the ID bytes
		  {
			int64_t n_tot_in_snap = 0;
			for(int k=0;k<6;k++)
			  n_tot_in_snap += header1.npart[k];

			if(dummy != (int)(ID_BYTES*n_tot_in_snap))
			  {
				fprintf(stderr,"ERROR: ID_BYTES =%zu and Ntot = %"PRId64" implies %"PRId64" bytes in the particle ids (in file `%s')  \n",
						ID_BYTES,n_tot_in_snap,ID_BYTES*n_tot_in_snap,snapshot_name);

				fprintf(stderr,"However I get bytes = %d with a corresponding ID_BYTES = %"PRId64" ..exiting\n",dummy,dummy/n_tot_in_snap);
				exit(EXIT_FAILURE);
			  }
		  }//end of scope

		  for(int k=0; k<6; k++) 
			fseek(fdmass, (long) (header1.npart[k]*ID_BYTES), SEEK_CUR);//casting to 'long' to avoid the signed/unsigned compilation warning.
		  my_fread(&dummy, sizeof(dummy), one, fdmass);
		  my_fread(&dummy, sizeof(dummy), one, fdmass);//set it read to read in the first of the particle masses
		}

	  
	  for(int k=0; k<6; k++)
		{ 
		  for(int npart=0; npart<header1.npart[k]; npart++)
			{ 
			  //read in the position from the snapshot
			  my_fread(&pos[0],sizeof(float),3*one,fd);

			  //compute histogram 3-d bin indices
			  xi = (int)(pos[0]/binsize);
			  yi = (int)(pos[1]/binsize);
			  zi = (int)(pos[2]/binsize);
			  if(xi >=0 && xi < Nbins && yi >=0 && yi < Nbins && zi>=0 && zi < Nbins)
				{
				  if(header1.mass[k] > 0.0)
					partmass = header1.mass[k];
				  else
					my_fread(&partmass,sizeof(float),one,fdmass);//read in the particle individual mass and add it to the histogram

				  histogram[xi][yi][zi] += partmass;
				  total_part_mass += partmass;
				}
			  else
				{
				  fprintf(stderr,"ERROR: (xi,yi,zi) = (%d,%d,%d) are not within [0,%d] for binsize = %lf \n",xi,yi,zi,Nbins,binsize);
				  fprintf(stderr,"pos[0] = %f pos[1] = %f pos[2] = %f ..exiting\n",pos[0],pos[1],pos[2]);
				  interrupted=1;//doesn't really matter -> exiting
				  exit(EXIT_FAILURE);
				}

			  index++;

			  //If you want to speed-up the code slightly -> put this entire chunk of 
			  //progress bar stuff outside of this loop. The progress bar will jump a lot
			  //but you will save quite a few computations. From my test on a 512^3, the
			  //code takes 35 seconds vs 38 seconds. Your call. 
			  if(PRINTSTEP > 0 && SMALLPRINTSTEP > 0)
				{ 
				  if(index%PRINTSTEP == 0)
					{ 
					  percent=(int)ceil((float)(index/PRINTSTEP))*10;
					  if(percent > 0)
						{ 
						  if(interrupted == 0)
							fprintf(stderr,"%02d%%",(int) percent);

						  my_snprintf(&(progressbarstring[end_of_string_index]),MAXLEN-end_of_string_index,"%02d%%",percent);
						  end_of_string_index += 3; //two digits for the number and one for the % symbol                                                                                                                
						}
					}
				  else
					{ 
					  if(index%SMALLPRINTSTEP==0)
						{ 
						  if(interrupted == 0)
							fprintf(stderr,".");
						  progressbarstring[end_of_string_index++] = '.';
						}
					}
				  progressbarstring[end_of_string_index+1] = '\0';
				  if(interrupted == 1)
					{ 
					  fprintf(stderr,"\n%s",progressbarstring);
					  interrupted = 0;
					}
				}

			}  
	  
		}
	  fprintf(stderr,"\nReading file `%s' ...done\n",snapshot_name);
	  interrupted=1;
	  fclose(fd);

	  if(FLAG_INDIVIDUAL_PARTMASS==1)
		fclose(fdmass);

	  if(nfiles > 1)
		{
		  t_sectionend = time(NULL);
		  //no need to check this string copy -> something else would have failed before this statement overflows
		  my_snprintf(outfname,MAXLEN,"Time to read in snapshot file# %d out of %d total files",ifile,nfiles);

		  print_time(t_sectionstart,t_sectionend,outfname);
		  interrupted=1;
		}
	  ifile++;
	}

  //now convert the mass in the bins into a normalized number density.
  for(int i=0;i<Nbins;i++)
    for(int j=0;j<Nbins;j++)
      for(int k=0;k<Nbins;k++)
	{
	  histogram[i][j][k] /= total_part_mass;//number count in the bins
	  histogram[i][j][k] /= (binsize*binsize*binsize);//co-moving number density in the bins
	}

  fprintf(stderr,"Calculating the histogram for Gadget snapshots ...done\n");
  //write out the histogram
  {
    my_snprintf(outfname,MAXLEN,"%s/histogram_numbers_%03d.f77binary",outpath,snapshot);
    fd = my_fopen(outfname,"w");
    int64_t N = Nbins*Nbins*Nbins;
    size_t blklen;
    
    fprintf(stderr,"\n\nWriting the histogram to `%s' ...",outfname);
    
#define BLKLEN my_fwrite(&blklen, sizeof(blklen), one, fd);
    //write out Nbins
    blklen = sizeof(int);
    BLKLEN;
    my_fwrite(&Nbins,sizeof(int),one,fd);
    BLKLEN;
    
    //write out the Gadget header
    blklen = sizeof(struct io_header);
    BLKLEN;
    my_fwrite(&header,sizeof(struct io_header),one,fd);
    BLKLEN;
    

    //write out the actual histogram data
    blklen=N*sizeof(double);
    BLKLEN;
    for(int i=0;i<Nbins;i++)
      for(int j=0;j<Nbins;j++)
	for(int k=0;k<Nbins;k++)
	  my_fwrite(&histogram[k][j][i],sizeof(histogram[k][j][i]),one,fd);
    
    BLKLEN;
    
    fclose(fd);
    fprintf(stderr,"..done\n");
    
    if(sizeof(size_t) != 4)
      fprintf(stderr,"To read in the histogram binary file in fortran, use `-frecord-marker=%zu' during compilation\n",
	      sizeof(size_t));//frecord-marker can only be 4 or 8. I could avoid this call to sizeof!
  }
  
  //free the histogram
  for(int i=0;i<Nbins;i++)
    {
      for(int j=0;j<Nbins;j++)
        free(histogram[i][j]);
      
      free(histogram[i]);
    }
  free(histogram);
  
  fprintf(stderr,"\n\n");
  t_codeend = time(NULL);
  print_time(t_codestart,t_codeend,"generate histogram");
  fprintf(stderr,"\n\n");
  return EXIT_SUCCESS;
}
