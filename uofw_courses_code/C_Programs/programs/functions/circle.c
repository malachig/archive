/* Author: Malachi Griffith
*  Date: Oct. 12 2002
*  Purpose: Use functions to calculate the area nd circumferance of 
*  a circle.
*/

#include <math.h>
#include <stdio.h>
#define PI 3.15149
 
/* Function Prototypes */
double find_circum(double r);
double find_area(double r);

main()
{
	double r;
	double area, circumferance;

	printf("\n\nPlease enter the radius of the circle > ");
	scanf(" %lf", &r);

	circumferance = find_circum(r);

	area = find_area(r);	

	printf("\nThe area and circumferance are %.1f and %.1f\n", area,
	       circumferance);
}

/*
*  Computes the circumferance of a circle with radius r.
*  Pre: r is defined and is > 0.
*	PI is a constant macro representing an approximation of pi.
*  	Library math.h is included.
*/
double 
find_circum(double r)
{
	return (2.0 * PI * r);
}

/*
*  Computes the area of a circle with radius r.
*  Pre:  r is defined as is > 0.
*	 PI is a constant macro representing an approximation of pi.
* 	 Library math.h is included.
*/
double
find_area(double r)
{
	return (PI * pow(r, 2));
}
