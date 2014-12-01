/* for.c */

/* Author: Malachi Griffith
*  Date: Sept. 20 2002 
*  Purpose: Counter-controlled repetition using the for structure.
*/

#include <stdio.h>

main()
{
	int counter;

	/* Initialisation, repetition condition and increment are 
	*  all included in the for structure header.
	*/

	for (counter = 0; counter < 10; counter++)
		printf(" %d\n", counter + 1);
}



