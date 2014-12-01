/* ifelse3.c */

/* Author: Malachi Griffith
*  Date: Sept. 28 2002 
*  Purpose: If else statements (Problems from pg 162 of text).
*/

#include <stdio.h>

main()
{
	/* Question 1-c  */
	int zero_count = 0;
	int minus_sum = 0;
	int plus_sum = 0;
	int number; 
        int count;

	for (count = 1; count <= 10; count++)
		{		
		printf("\nPlease enter an integer number > ");
		scanf(" %d", &number);
		if (number == 0)
			zero_count++; 
		if (number > 0)
			plus_sum += number;	
		if (number < 0)
			minus_sum += number;
		}	
	printf("\n\nThe number of 0 entries was: %d", zero_count);
	printf("\nThe sum of the positive numbers was: %d", plus_sum);
	printf("\nThe sum of the negative numbers was: %d\n", minus_sum);
	
}



