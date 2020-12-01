/* score2.c */

/* Author: Malachi Griffith
*  Date: Oct. 6 2002
*  Purpose: Compute the sum of a list of exam scores 
*	USE OF A FOR STATEMENT TO IMPLEMENT A SENTINLE CONTROLLED LOOP.
*/

#include <stdio.h>
#define SENTINEL -99

int
main(void)
{
	int sum = 0,	/* Output - sum of scores input so far */
	    score;	/* Input - current score */

	/* Accumulate sum of all scores */
	printf("Enter first score (or %d to quit) > ", SENTINEL);
	scanf("%d", &score);

	for (scanf("%d", &score);
	     score != SENTINEL;
	     scanf("%d", &score))	
	{
	sum += score;
	printf("Enter next score (%d to quit) > ", SENTINEL);
	}

	printf("\nSum of exam scores is %d\n", sum);

	return(0);
}



