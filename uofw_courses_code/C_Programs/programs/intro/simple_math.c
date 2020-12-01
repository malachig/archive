/* simple_math*/
/* Author: Malachi Griffith
*  Date: Sept. 20 2002
*  Purpose:
*/

#include <stdio.h>

main()
{
	/* Program Inputs */
	double 	x, y;

	/* Program Outputs */
	double	sum,
		difference,
		product,
		quotient;
	printf("\nType in the first integer > ");
	scanf(" %lf", &x);
	printf("Type in the second integer > ");
	scanf(" %lf", &y);

	sum = (double) (x + y);
	printf("The sum of x and y is: %5.2f\n", sum);
	difference = (double) (x - y);
	printf("The difference between x and y is: %5.2f\n", difference);
	product = (double) (x * y);		
	printf("The product of x and y is: %5.2f\n", product);
	quotient = (double) (x / y);
	printf("The quotient of x and y is: %5.2f\n", quotient);
}
