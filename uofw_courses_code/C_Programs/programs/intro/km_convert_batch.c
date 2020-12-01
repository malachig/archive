/*km_convert_batch.c*/
/*
* Converts distances from miles to kilometers
*/

#include <stdio.h>  /* printf, scanf definitions */
#define kms_per_mile 1.609  /* conversion constant */

int
main(void)
{
	double miles, /* diustance in miles*/
		kms; /* equivalent distance in kilometers */ 

	/* Get and echo the distance in miles */
	scanf ("%lf", &miles);
	printf("The distance in miles is %.2f.\n", miles);

	/*convert the distance to kilometers. */
	kms = kms_per_mile * miles;
 
	/*display the distance in kilometers. */
	printf("That equals %f kilometers. \n", kms);

	return (0);
}





