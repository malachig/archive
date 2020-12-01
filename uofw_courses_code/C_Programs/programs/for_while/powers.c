/* powers.c */

/* Author: Malachi Griffith
*  Date: Oct. 5 2002 
*  Purpose:  Illustrate the use of a simple while loop
*/

#include <stdio.h>
#include <math.h>

main()
{
	int counter = 0;
	int power;

	while (counter <= 6)
	{
	power = (pow(2, counter));
	printf("\n\t %d\t %d", counter, power);
	counter++;
	}
	printf("\n\n");
}



