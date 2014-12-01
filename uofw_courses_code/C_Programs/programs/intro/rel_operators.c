/* rel_operators.c*/
/* Author: Malachi Griffith 
*  Date: Sept. 17 2002
*  Purpose: illustrate the use of relational operators.
*/

/* The relational operators are defined as follows:
*	> means greater than,
*	>= means greater than or equal to,
*	< means less than,
*	<= means less than or equal to,
*  These four all have the same precedence.  The following two have
*  less precedence:
*	== means equal to,
*	!= means not equal to.
*/

#include <stdio.h>

main()
{
	int i = 10;
	int j = 10;
	int k = 15;

	printf("i==j yields: %d\n", i==j);
	printf("i!=j yields: %d\n", i!=j);
	printf("i==k yields: %d\n", i==k);
	printf("i!=k yields: %d\n", i!=k);
	printf("i<j yields: %d\n", i<j);
	printf("i<k yields: %d\n", i<k);

/* Note: The output for each of these statements will be 1 or 0,
*  where 1 means TRUE and 0 means FALSE
*/

}


