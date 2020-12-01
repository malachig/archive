/* Author: Malachi Griffith
*  Date: Nov. 2 2002 
*  Purpose: Illustrate use of indirection (*) and address of operators (&)
*/

#include <stdio.h>

main()
{
	int i;
	int *i_ptr = &i;

	printf("\nThe address of i is %p", &i);
	printf("\nThe value of i_ptr is %p", i_ptr);

	/* Place a value into i using indirection */

	*i_ptr = 17;

	printf("\n\nThe value of i is   %d", i);
	printf("\nThe value of *i_ptr is %d", *i_ptr);
}



