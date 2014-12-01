/* Author: Malachi Griffith
*  Date: Nov. 10 2002 
*  Purpose: Displays the values along the diagonal of a 10 by 10 matrix.
*/

#include <stdio.h>
#include <math.h>

main()
{
	int value;
	int matrix[10][10];
	int row, col;
	int n = 9;
	printf("\n");

	for (row = 0; row < 10; row++)
		for (col = 0; col < 10; col++)
			matrix[row][col] = rand();

	for(row = 0; row < 10; row++)
	{
		for (col = 0; col < 10; col++)
			printf("%d\t", (matrix[row] [col])/100);

	printf("\n\n");
	}
	
	for(row = 0; row < 10; row++)
	{
		for (col = 0; col < 10; col++)
		{
			if(col == row)
				printf("%d", (matrix[row] [col])/100);
			else if (col == n)
				printf("\t%d", (matrix[row] [col])/100);
			else
				printf("\t");	
		}
	n--;	
	printf("\n");
	}
	




}



