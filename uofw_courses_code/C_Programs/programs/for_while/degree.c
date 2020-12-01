/* degree.c */

/* Author: Malachi Griffith
*  Date: Oct. 5 2002
*  Purpose: Displays a table of celsius to farenheit conversions
*           Uses a loop that counts by units of 5 or whatever I set
*	    CSTEP to.
*/

#include <stdio.h>

#define CBEGIN 10
#define CLIMIT -5
#define CSTEP 5

int
main()
{
	/* Variable declarations */
	int celsius;
	double farenheit;

	/* Display the table heading */
	printf("   Celsius    Fahrenheit\n");

	/*Display the table */
	for (celsius = CBEGIN;
	     celsius >= CLIMIT;
	     celsius -= CSTEP)
	{
	farenheit = 1.8 * celsius + 32.0;
	printf("%6c%3d%8c%7.2f\n", ' ', celsius, ' ', farenheit);
	}
	
	return(0);
}



