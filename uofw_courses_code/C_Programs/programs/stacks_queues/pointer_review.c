/* Author: Malachi Griffith
*  Date: Dec. 7 2002 
*  Purpose: Function with pointers as output parameters.
*/

#include <stdio.h>

void long_division(int dividend, int divisor, int *quotientp,
		   int *remainderp);

int
main(void)
{
	int quot, rem;

	long_division(40, 3, &quot, &rem);
	printf("40 divided by 3 yields quotient %d ", quot);
	printf("and remainder %d\n", rem);
	
	return(0);
}

/*
*  Performs long division of two integers, storing quotient 
*  in variable pointed to by quotientp and remainder in 
*  variable pointed to by remainderp.
*/
void
long_division (int dividend, int divisor, int *quotientp,
	       int *remainderp)
{
	*quotientp = dividend / divisor;
	*remainderp = dividend % divisor;
}
