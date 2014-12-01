/* Author: Malachi Griffith
*  Date: Nov. 2 2002
*  Purpose: Illustrate difference between simple arguments and pointers. 
*  Using Functions but not really working.
*/

#include <stdio.h>

void exchange (int x, int y)
{
	int temp;
	temp = x;
	x = y;
	y = temp;

	printf("in function x is %d, y is %d\n", x, y);
}

main()
{
	int x = 3;
	int y = 4;
	int temp;

	printf("in main, x is %d, y is %d\n", x, y);

	exchange(x, y);

	printf("in main, x is %d, y is %d\n", x, y);
}
