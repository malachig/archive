/*km_convert_pointer.c*/
/* Author: Malachi Griffith
*  Date:  Sept. 21 2002
*  Purpose: To illustrate the use of program-controlled input and output
*	    files.  Use of file pointer variables and the functions fopen,
*	    fscanf, and fprintf will be illustrated
*/

#include <stdio.h>
#define kms_per_mile 1.609 /* Conversion constant */

int 
main(void)
{
	double miles,
	       kms;

	FILE   *inp,     /* Pointer to input file */
	       *outp;    /* Pointer to output file */

	/* Open the input and output files. */
	inp = fopen("km_convert.dat", "r");
	outp = fopen("km_convert.out", "w");

	/* Get and echo the distance in kilometers */
	fscanf(inp, "%lf", &miles);
	fprintf(outp, "The distance in miles is %.2f. \n", miles);

	/* Convert the distance to kilometers. */
	kms = kms_per_mile * miles;

	/* Display the distance in kilometers. */
	fprintf(outp, "That equals %.2f kilometers. \n", kms);

	/* Close the files */
	fclose(inp);
	fclose(outp);

	return(0);
/*  When this program is compiled and a.out is run, it reads the file,
*   km_convert.dat for input and generates its output to the file,
*   km_convert.out.
*/

}

