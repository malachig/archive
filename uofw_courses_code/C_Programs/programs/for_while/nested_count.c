/* nested_count.c */

/* Author: Malachi Griffith
*  Date:  Oct. 5 2002
*  Purpose: Illustrates a pair of nested counting loops.
*/

#include <stdio.h>

int
main(void)
{
	int i, j;	/*loop control variables */
	
	printf("	I	J\n");	/* prints column labels */

	for (i = 1; i < 4; ++i)
	{
		printf("Outer %3d\n", i);
		
		for (j = 0; j < i; ++j)
		{
		printf("  Inner%10d\n", j);
		}
	}
	return(0);
}



