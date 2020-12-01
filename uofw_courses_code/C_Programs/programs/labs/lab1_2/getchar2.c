/* Author: Malachi Griffith
*  Date: Sept. 16 2002
*  Purpose: Illustrate use of the "getchar" command
*/

#include <stdio.h>
main()
{
	char c1;
	char c2;
	
	printf("Enter a 2-letter initials \n");

	c1 = getchar();
	c2 = getchar();

	printf("Your initials are %c and %c >\n", c2, c1);
}

