/* dowhile_input.c */

/* Author: Malachi Griffith
*  Date:  Oct.6 2002
*  Purpose: Input loop that scans pair of integers until it finds a
*  pair where the first integer evenly divides the second.
*/

#include <stdio.h>

main()
{
	int first,
	    second;
	int test;

	do
	{
	printf("\n\nEnter the first integer > ");	
	scanf("%d", &first);
	printf("Enter the second integer > " );
	scanf("%d", &second);
	test = second % first;
	}
	while (test != 0);
	printf("\nThe first evenly divides the second\n\n");
}



