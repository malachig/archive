/* average.c */

/* Author: Malachi Griffith
*  Date:  Sept. 28 2002 
*  Purpose:  Uses a for loop to calculate the average of a set of n numbers.
*/

#include <stdio.h>

main()
{
	int n;
	int value;
	int total = 0;
	int average;
	int counter;

	printf("\n\nEnter the number of values you wish to enter > ");	
	scanf(" %d", &n);

	for (counter = 0; counter < n; counter++)
	{
	printf("Enter an integer value > ");
	scanf(" %d", &value); 
	total += value;
	}
	average = total / n;
	printf("\n The average of the entered values is: %d\n\n", average);
}



