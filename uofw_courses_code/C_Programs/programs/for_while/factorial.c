/* factorial.c */

/* Author: Malachi Griffith
*  Date: Oct. 5 2002
*  Purpose: Computes n a factorial
*/

#include <stdio.h>

int
factorial(int n)
{
	int i,		/* local variables */
	product;	/* accumulator for product computation */

	product = 1;
	/* Computes the product n x (n-1) x (n-2) x ... x 2 x 1 */

	for (i = n; i > 1; --i)
	{
		product = product * i;
	}

	/*returns function result */
	return (product);
}

main()
{
	int x;
	int answer;
	printf("\n\nplease enter the factorial you wish to find > ");
	scanf("%d", &x);

	answer = factorial(x);
	printf("\nThe factorial of that number is : %d\n\n", answer);
}

