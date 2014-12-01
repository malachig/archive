/* break.c */

/* Author: Malachi Griffith
*  Date: Sept. 20 2002
*  Purpose: To illustrate the use of the break command
*  The break command causes immediate exit from any loop structure.
*/

#include <stdio.h>
#define START 1
#define END 10
#define EXCEPTION 5

main()
{
	int x;

	for(x = START; x <= END; x++)
	{
		if(x == EXCEPTION)
			break;

		printf("%d ", x);
	}

	printf("\nBroke out of loop at x == %d\n", x);
}



