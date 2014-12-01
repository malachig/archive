/* Author: Malachi Griffith
*  Date: Nov. 2 2002
*  Purpose: Illustrate difference between simple arguments and pointers. 
*  NOT using Functions.
*/

#include <stdio.h>

main()
{
	int x = 3;
	int y = 4;
	int temp;

	printf("x is %d, y is %d\n", x, y);
	temp = x;
	x = y;
	y = temp;
	printf("x is %d, y is %d\n", x, y);
}
