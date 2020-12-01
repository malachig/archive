/* Author: Malachi Griffith
*  Date: Nov. 2 2002 
*  Purpose: Illustrate the use of pointer arithmetic to print a string.
*/

#include <stdio.h>

main()
{
	char *s = "Here is a string";

	while(*s != '\0')  /* The '\0' is the end of string mark */
	{
	putchar(*s);
	s++;  		/* increment the address of the pointer */
	}

	printf("\n");
}
