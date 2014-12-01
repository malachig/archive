/* char_assignments.c */

/* Author: Malachi Griffith
*  Date: Sept. 20 2002
*  Purpose: To illustrate the use of assignment operators and logical statements when working
*	    with characters instead of numbers.
*/

#include <stdio.h>

main()
{

	int 	
		uppercase,
		divisor,
		m,
		n;
	char	ch;

	printf("\n\nPlease enter a single letter > ");
	scanf(" %c", &ch);
	uppercase = (ch >= 'A' && 'Z' >= ch);
	printf("Test of whether this character is an uppercase letter: %d\n\n", uppercase);

	printf("\n\nPlease enter an integer value for m and n > ");
	scanf(" %d%d", &m, &n);
	divisor = ((n % m) == 0);
	printf("Test of whether m is a divisor of n: %d\n\n", divisor);

}
