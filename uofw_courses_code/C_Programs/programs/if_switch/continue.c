/* continue.c */

/* Author: Malachi Griffith
*  Date: Sept. 20 2002 
*  Purpose: A program to demonstrate the use of continue.
*  In a while, for, or do\while statement the continue statement
*  skips the remaining statements in the body of the structure and performs 
*  the next iteration of the loop.
*/

#include <stdio.h>
#define START 1
#define END 10
#define EXCEPTION 5
main()
{
	int x;

	for (x = START; x <= END; x++)
	{
		if(x == EXCEPTION)
			continue;

		printf("%d ", x);
	}

	printf("\nUsed continue to skip printing the value 5\n\n");
}



