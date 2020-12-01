/* Author: Malachi Griffith
*  Date: Oct. 7 2002
*  Purpose: A program for testing the function scale().
*/

#include <stdio.h>
#include <math.h>

double scale(double x, int n);

main()
{
	int num_2;
	double num_1;
	
	printf("Enter two integers, a number and then x where 10eX");	
	scanf("%lf %d", &num_1, &num_2);
	printf("Result of call to function scale is %.2f\n",
	       scale(num_1, num_2));
}

/* A function that multiplies its first argument by the power
*  of 10 specified by its second argument.  It uses the standard
*  math library function pow().
*/

double scale(double x, int n)
{
	double scale_factor;
	scale_factor = pow(10, n);

	return (x * scale_factor);
}



