/* Author: Malachi Griffith
*  Date: Oct. 7 2002
*  Purpose: A program for finding the largest of three integers. 
*/

#include <stdio.h>

/* Function prototype */
int maximum (int, int, int);

main()
{
	int a;
	int b;
	int c;

	printf("Please enter three integers: ");
	scanf("%d%d%d", &a, &b, &c);
	printf("\nMaximum is: %d\n", maximum (a, b, c));
}

int
maximum(int x, int y, int z)
{
	int max = x;

	if (y > max)
		max = y;

	if (z > max)
		max = z;

	return(max);
}



