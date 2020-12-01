/* Author: Malachi Griffith
*  Date: Oct. 7 2002 
*  Purpose: A program that uses a user defined function to print 
*	    a hollow square of stars.
*/

#include <stdio.h>

/* Function prototype */
void square(int);

main()
{
	int size;
	
	printf("Please enter size of square: ");
	scanf("%d", &size);
	
	square(size);
}

void square(int size)
{
	int row;
	int col;
	
	printf("\n");
	
	for(col = 0; col < size; col++)
		printf("*");
	printf("\n");

	for (row = 0; row < size - 2; row++)
	{
		printf("*");
		for (col = 0; col < size - 2; col++)
			printf(" ");
		printf("*\n");
	}
	if (size > 1)
		for(col = 0; col < size; col++)
			printf("*");
	printf("\n");
	return;
}



