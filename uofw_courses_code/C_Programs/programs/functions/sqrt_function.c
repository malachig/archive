/* Author: Malachi Griffith
*  Date:  Oct. 8 2002
*  Purpose: Performs three square root computations
*/

#include <stdio.h>
#include <math.h>

int
main(void)
{
	double first, second,
	       first_sqrt,
	       second_sqrt,
	       sum_sqrt;

	/* Get first number and display its square root */
	printf("Enter the first number > ");
	scanf("%lf", &first);
	first_sqrt = sqrt(first);
	printf("The square root of the first number is %.2f\n", first_sqrt);

	/* Get second number and display its square root. */
	printf("Enter the second number > ");
	scanf("%lf", &second);
	second_sqrt = sqrt(second);
	printf("The square root of the second number is %.2f\n", second_sqrt);

	/* Display the square root of the sum of the two numbers */
	sum_sqrt = sqrt(first + second);
	printf("The square root of the sum of the two numbers is %.2f\n", sum_sqrt);

	return(0);
}
