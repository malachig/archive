/* Author: Malachi Griffith
*  Date: Oct. 7 2002
*  Purpose: A program for calculating the square of a number 
*/

#include <stdio.h>

int square(int a)
{
	int b;
	b = a * a;
	return (b);
}

main()
{
	int x;
	
	printf("Please enter an integer: ");
	scanf("%d", &x);
	printf("The square of %d is %d\n", x, square(x));
	printf("The square of 3 times %d is %d\n", x, square(3 * x));
}
