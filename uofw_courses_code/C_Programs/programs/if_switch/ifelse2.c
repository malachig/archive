/* ifelse2.c */

/* Author: Malachi Griffith
*  Date: Sept. 28 2002 
*  Purpose: If else statements (Problems from pg 162 of text).
*/

#include <stdio.h>

main()
{
	int product = 5;
	int item;
	int abs_dif,
	    x,
	    y;

	/* Question 1-a */	
	printf("\n\nPlease enter a integer value > ");
	scanf(" %d", &item);

	if (item != 0)
		printf("\n\nThe product of 5 and this number is: %d", product * item);
	else
		printf("\nMultiplication by 0");

	/* Question 1-b */
	printf("\n\nPlease enter a value for x > ");
	scanf(" %d", &x);
	printf("Please enter a value for y > ");
	scanf(" %d", &y);

	if (x >= y)
		abs_dif = x - y;
	else
		abs_dif = y - x;
	printf("\nThe absolute difference between these two values is: %d\n\n", abs_dif);

	/* Question 1-c 
	int zero_count = 0;
	int minus_sum = 0;
	int plus_sum = 0;
	int number; 
        int count;

	for (count = 1; count <= 10; count++)
		
		printf("\nPlease enter an integer number");
		scanf(" %d", &number);
		if (number = 0)
			zero_count++; 
		if (number > 0)
			plus_sum += number;	
		if (number < 0)
			minus_sum += number;
		
	printf("\n\nThe number of 0 entries was: %d", zero_count);
	printf("\nThe sum of the positive numbers was: %d", plus_sum);
	printf("\nThe sum of the negative numbers was: %d", minus_sum);
	*/
}



