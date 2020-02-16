// *****************************************************************
//  
//    Author:           Ehsan Kourkchi
//  
//    DATE:             Jan 3 2010
//  
//    FILE:             pixval.cpp     
//  
//    DESCRIPTION:      Given a fits file, this reads the pixel value at the given x,y coordinate
//  
//  ****************************************************************/
// Usage:
//   $ pixval  <imgae.fits>  x y
//    where <image.fits> is a 2-dimentional fits file
//    x and y are the coordinates at which the pixel value is desired.
//  example: pixval  segmentation.fits  50 70
// 
//   To compile this program use:
//   $ c++ -o pixval pixval.cpp -I<where cfitsio installed>/cfitsio -Wno-deprecated -lcfitsio
//
//  ****************************************************************/


#include <string.h>
#include <stdio.h>
#include "fitsio.h"
#define pi 4*atan(1.)
#include "stdlib.h"
#include "time.h"
#include "math.h"


int main(int argc, char *argv[])
{
    fitsfile *fptr,*outfptr;   /* FITS file pointer, defined in fitsio.h */
    int status = 0;   /* CFITSIO status value MUST be initialized to zero! */
    int bitpix, naxis, ii, anynul, jj;
    long naxes[2] = {1,1}, fpixel[2] = {1,1};
    double *pixels,*outpix;
    char format[20], hdformat[20];
    double x0,y0,A,B,theta,r_kron,no,A0,B0;
    double x,y,X_m,Y_m,check;
    double coef;
    if (argc != 4) {
      printf("\n\n\n\n\nError: Please enter the arguments correctly ... \n");
      printf("Example: $ pixval  <imgae.fits>  x y ... \n");
      printf("where <image.fits> is a 2-dimentional fits file \n");
      printf("wx and y are the coordinates at which the pixel value is desired. \n");
      printf("Ehsan Kourkchi 3 Jan 2010 \n \n");
      return(0);
    }

ii=atoi(argv[2]);
fpixel[1]=long (atoi(argv[3]));



    if (!fits_open_file(&fptr, argv[1], READONLY, &status))
    {
        if (!fits_get_img_param(fptr, 2, &bitpix, &naxis, naxes, &status) )
        {
          if (naxis > 2 || naxis == 0 || naxis == 1)
             printf("Error: only 2D images are supported\nMulti extension images are not supported.\n");
          else
          {


	if (ii<naxes[0] && ii>0 &&  fpixel[1]>0 &&  fpixel[1]<naxes[1]){ 
            pixels = (double *) malloc(naxes[0] * sizeof(double));

            if (pixels == NULL) {
                printf("Memory allocation error\n");
                return(1);
            }


               fits_read_pix(fptr, TDOUBLE, fpixel, naxes[0], NULL,
                    pixels, NULL, &status);  /* read row of pixels */
                 

               
printf("%lf\n",pixels[ii]); // The desired pixel value
              

            free(pixels);
            

          } else printf("0.000\n"); }
        }
        fits_close_file(fptr, &status);
    } 

    if (status) fits_report_error(stderr, status); /* print any error message */
    return(status);
}
