/* multiply.c */

/* Author: Malachi Griffith
*  Date: Oct. 5 2002
*  Purpose: Multiply a list of numbers until product exceeds 10000
*/

#include <stdio.h>

main()
{
	double product = 1;
	double number;
	
	while (product < 10000)
	{
		printf("\n%.0f\n", product);
		printf("Enter a number to multiply > ");
		scanf("%lf", &number);
		product = product * number;
	}
}



