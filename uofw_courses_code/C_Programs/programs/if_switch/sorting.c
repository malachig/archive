/* sorting.c */

/* Author: Malachi Griffith
*  Date: Sept. 28 2002 
*  Purpose: Sorts 2 numbers so that the smaller is stored in x and 
*           the larger is stored in y even if they start out the other way.
*/

#include <stdio.h>

main()
{
	int x, y, temp;

	printf("\n\nEnter a value for x > ");
	scanf(" %d", &x);	
	printf("Enter a value for y > ");
	scanf(" %d", &y);

	if (x > y)
		{
		temp = x;
		x = y;
		y = temp;
		}
	printf("\nThe sorted values of x and y are: %d and %d\n\n", x, y);
}



