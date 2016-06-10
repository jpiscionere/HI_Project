/* headr -- read header info from a fastfood.out file
   headr
      < fastfood.out file
      > header info

      -- like DHW's header, but aslo reads the 1st zstep 
*/
#include <stdio.h>
#include <stdlib.h>

int main(int argc,char **argv)
{
    int idat[5] ;
    float xdat[9] ;
    float znow;
    
    ftread(idat,4,5,stdin) ;
    ftread(xdat,4,9,stdin) ;
    ftread(&znow,4,1,stdin);

    printf("n = %d, np = %d , iskip = %d\n",idat[0],idat[1],idat[2]) ;
    printf("rcube = %f , h = %f\n",xdat[0],xdat[1]) ;
    printf("xindx = %f, rsm = %f, itrans = %d\n",xdat[2],xdat[3],idat[3]) ;
    printf("rgaus = %f, rth = %f, delnorm = %f\n",xdat[4],xdat[5],xdat[6]) ;
    printf("gnu = %f, eta = %f, jseed = %d\n",xdat[7],xdat[8],idat[4]) ;
    printf("znow = %f\n",znow);
    exit(EXIT_SUCCESS);
}
