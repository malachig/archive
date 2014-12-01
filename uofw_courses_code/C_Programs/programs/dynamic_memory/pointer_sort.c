/* Author: Malachi Griffith
*  Date: Nov. 21 2002 
*  Purpose: Illustrate simple exchanging of pointers
*/

#include <stdio.h>
#include <stdlib.h>  /* Needed for dynamic memory allocation */

typedef int *intptr;

void exchange(intptr *x, intptr *y)
{
	intptr temp;
	temp = *x;
	*x = *y;
	*y = temp;
}

main()
{
	int x = 3, y = 4;
	
	intptr p1, p2;
	
	p1 = &x;
	p2 = &y;

	printf("Before exchanging ptrs, *p1 = %d, *p2 = %d\n", *p1, *p2);

	exchange(&p1, &p2);

	printf("After exchanging ptrs, *p1 = %d, *p2 = %d\n", *p1, *p2);
}
