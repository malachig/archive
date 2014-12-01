/* Author: Malachi Griffith
*  Date:  Oct. 6 2002 
*  Purpose: Displays a multiplication table for numbers 0 to 9
*/

#include <stdio.h>

int
main(void)
{
	int column,
	    row;
	printf("\n\n");
	
	for (column = 0; column <= 4; column++)
	{
		printf("\n");	
		for (row = 0; row <= column; row++)
		{
			printf("%2d", row);
		}
	}
	
	for (column = 3; column >= 0; column--)
	{
		printf("\n");
		for (row = 0; row <= column; row++)
		{
			printf("%2d", row);
		}
	}

printf("\n\n");
	return(0);
}



