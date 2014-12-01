/* inch_cent.c */

/* Author: Malachi Griffith
*  Date: Oct. 5 2002
*  Purpose: Displays a table of inches to centimeter conversions.
*/

#include <stdio.h>

#define inch_to_centimeter 2.54

main()
{
	int inch_start;
	int inch_end;
	int increment;

	int value;
	double centimeters;
	
	printf("\n\nPlease enter the begining of the range > ");
	scanf("%d", &inch_start);
	printf("Please enter the end of the range > ");
	scanf("%d", &inch_end);
	printf("Please enter the increment you wish > ");
	scanf("%d", &increment);

	printf("\n\n\tInches\tCentimeters");

	for (value = inch_start;
	     value <= inch_end;
	     value += increment)
	{
	centimeters = (double) (value * inch_to_centimeter);
	printf("\n\t%d\t%.2f", value, centimeters);
	}
	printf("\n\n");
}
