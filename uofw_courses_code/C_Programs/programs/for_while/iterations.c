/* iterations */

/* Author: Malachi Griffith
*  Date:  Oct. 5 2002 
*  Purpose: To compare the addition of 1 + 2 + ... + n to the 
*  multiplication of (n * (n + 1)) / 2.
*/

#include <stdio.h>

main()
{
	double value1;
	double value2;	
	int iterations;
	int counter = 0;

	printf("\n\nPlease enter the number of iterations > ");
	scanf("%d", &iterations);

	while (counter <= iterations)
	{
	value1 = value1 + counter;
	counter++;
	}
	
	printf("\nThe value of (1 + 2 + ... + n) is %.0f\n", value1);
	value2 = ((iterations * (iterations + 1)) / 2);
	printf("\nThe value of ((n * (n + 1)) / 2) is %.0f\n", value2); 

	if (value1 == value2)
		printf("\n\nThe values are the same!\n");
	else
		printf("\n\nThe values are different!\n");
}



