/* nestedloop.c */

/* Author: Malachi Griffith
*  Date: Sept. 20 2002 
*  Purpose: Program which uses nested loops to calculate a multiplication table
*/

#include <stdio.h>
#define MAXROW 5
#define MAXCOL 10

main()
{

	int row;
	int col;

	printf("\n\n\t\t\tMultiplication Table");
	
	for(row = 1; row <= MAXROW; row++)
	{
		printf("\n\n");
		for(col = 1; col <= MAXCOL; col++)
			printf("%6d", row * col);
	}

	printf("\n\n");
}



