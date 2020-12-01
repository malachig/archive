/*km_convert.c*/
/*
* Converts distances from miles to kilometers
* Since this is a batch program you must redirect the input when you run 
* the program.  Compile by saying %cc km_convert_batch.c  then when you
* run the program, type % a.out <km_input.txt  
*/

#include <stdio.h>  /* printf, scanf definitions */
#define kms_per_mile 1.609  /* conversion constant */

int
main(void)
{
	double miles, /* diustance in miles*/
		kms; /* equivalent distance in kilometers */ 

	/* Get the distance in miles */
	printf ("Enter the distance in miles> ");
	scanf ("%lf", &miles);

	/*convert the distance to kilometers. */
	kms = kms_per_mile * miles;
 
	/*display the distance in kilometers. */
	printf("That equals %f kilometers. \n", kms);

	return (0);
}





