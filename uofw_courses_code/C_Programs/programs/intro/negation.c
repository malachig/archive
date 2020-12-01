/*negation.c*/
/*
*  Author: Malachi Griffith
*  Date: Sept. 17 2002
*  Purpose: illustrate the use of the logical operator, NEGATION.
*  The symbol "!" is used to specify negation.
*  It reverses the truth value of its operand.
*/

#include <stdio.h>

main()
{
	int k = 15;

	printf("Reversing 0 gives: %d\n", !0);
	printf("Reversing k gives: %d\n", !k);
	printf("Reversing (k == 15) gives: %d\n", !(k == 15));
	printf("Reversing (k == 10) gives: %d\n", !(k == 10));
}


