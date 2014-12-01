/* Author: Malachi Griffith
*  Date: Nov. 3 2002 
*  Purpose: Test function "order by ordering three numbers.
*  In this example the same parameter is used both to send a value to 
*  a function but also to return a value to main.
*/

#include <stdio.h>

/* Function Prototypes */
void order(double *smp, double *lgp);

int
main(void)
{
	double num1, num2, num3;  /* threed numbers to put in order */

	/* Gets test data */
	printf("Enter three numbers seperated by blanks > ");

	scanf("%lf%lf%lf", &num1, &num2, &num3);
	/* Orders the three numbers */
	order(&num1, &num2);
	order(&num1, &num3);
	order(&num2, &num3);

	/* displays results */
	printf("The numbers in ascending order are: %.2f %.2f %.2f\n",
		num1, num2, num3);
	
	return(0);
}

/* function order */
/* Arranges arguments in ascending order.
*  Pre: smp and lgp are addresses of defined type double variables
*  Post: variable pointed to by smp contains the smaller of the type 
*	 double values; varialbe pointed to by lgp contains the larger.
*/

void
order(double *smp, double *lgp)
{
	double temp;  /* temporary variable to hold one number during swap. */
	
	/* compares values pointed to by smp and lgp and swaps if neccessary */

	if (*smp > *lgp)
	{
		temp = *smp;
		*smp = *lgp;
		*lgp = temp;
	}
}
	

