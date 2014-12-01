/* Author: Malachi Griffith
*  Date: Oct. 24 2002 
*  Purpose: Calculates the ratio between successive speeds of a six-speed
*  gearbox.  Use formula 5th-root of (M/m) where M = max RPM, m = min RPM.
*/

#include <stdio.h>
#include <math.h>

/* Function Prototypes */
double speeds_ratio(int max_rpm, int min_rpm);

main()
{
	int max, min;
	double ratio;

	printf("\nEnter the maximum RPM > ");
	scanf("%d", &max);
	printf("Enter the minimum RPM > ");
	scanf("%d", &min);

	ratio = speeds_ratio(max, min);

	printf("\nThe ratio between successive speeds of a six speed ");
	printf("gearbox with maximum speed of %d rpm ", max);
	printf("and minimum speed of %d rpm is %lf.\n\n", min, ratio);	
}

/*
*  Function speeds_ratio
*/
double
speeds_ratio(int max_rpm, int min_rpm)
{
	double gear_ratio;
	double quotient;
	quotient = (double) (max_rpm/min_rpm);

	gear_ratio = pow(quotient, 0.2);
	return(gear_ratio);
}


