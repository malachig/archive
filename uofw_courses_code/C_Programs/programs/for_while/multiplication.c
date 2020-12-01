/* multiplication.c */

/* Author: Malachi Griffith
*  Date:  Oct. 6 2002 
*  Purpose: Displays a multiplication table for numbers 0 to 9
*/

#include <stdio.h>
#define MAX_COL 9
#define MAX_ROW 9

int
main(void)
{
	int column,
	    row;
	printf("\n\n");

	for (column = 0; column <= MAX_COL; column++)
	{
		printf("\n");	
		for (row = 0; row <= MAX_ROW; row++)
		{
			printf("%4d", row * column);
		}
	}
	printf("\n\n");
	return(0);
}



