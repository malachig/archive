/* gcd.c */

/* Author: Malachi Griffith
*  Date: Oct. 25 2002
*  Purpose: Finds the greatest common divisor of two integers.
*/

#include <stdio.h>
#include <math.h>

main()
{
	int num1, num2;
	int abs_num1, abs_num2; 
	int remainder;
	int divisor;
	int temp;
	int value;

	printf("\nEnter two integer values seperated by space > ", num1, num2);
	scanf("%d%d", &num1, &num2);

	abs_num1 = abs(num1);
	abs_num2 = abs(num2);
	
	divisor = abs_num2;
	value = abs_num1;
	remainder = value % divisor; 

	while (remainder != 0)
	{
	temp = divisor;	
	divisor = remainder;	
	value = temp;
	
	remainder = value % divisor;
	}
	printf("\nThe divisor is %d, value is %d, remainder is %d\n", 
		divisor, value, remainder);
	printf("\nTherefore the greatest common divisor is %d\n\n", divisor);
}



