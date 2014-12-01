/* digits.c */

/* Author: Malachi Griffith
*  Date: 
*  Purpose:
*/

#include <stdio.h>

main()
{
	int number;
	int digit;
	int value;
	int i;
	
	printf("\nEnter an integer value > ");
	scanf("%d", &number);

	for (i = 1; i <=10; i++)	
	{
		value = number / 10;
		digit = number % 10;
		number = value;
		if (value < 10 && digit == 0)
			break;
		printf("\n%d", digit);
	}

	printf("\n\n");





}



