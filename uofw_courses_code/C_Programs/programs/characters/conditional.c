/* Author: Malachi Griffith
*  Date: Oct. 14 2002
*  Purpose: A program for demonstrating characters in conditional statements
*/

#include <stdio.h>
#define TEST 'q'

main()
{
	char ch;
	int i;

	for (i = 0; ch != TEST; i++)
	{
		ch = getchar();
		putchar(ch);
	}

	printf("\nThe character %c is found in position %d.\n", ch , i);
}
