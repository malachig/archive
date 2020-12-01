/* Author: Malachi Griffith
*  Date: Dec. 7 2002  
*  Purpose: Simply to illustrate some examples of how
*  to use the malloc() function to dynamically create
*  memory for particular data pointer types.
*/

#include <stdio.h>
#include <stdlib.h>

main()
{
	/* Without Dynamic Allocation */
	int *nump;
	char *letp;
	planet_t *planetp;

	/* With Dynamic Allocation */

	nump = (int *)malloc(sizeof(int));
	letp = (char *)malloc(sizeof(char));
	planetp = (planet_t *)malloc(sizeof(planet_t));


	/* Storing values in the allocated memory - use indirection
	   operator. */

	*nump = 307;
	*letp = 'Q';
	*planetp = blank_plantet;	

}



