/* Author: Malachi Griffith
*  Date: Oct. 14 2002
*  Purpose: A program for demonstrating the use of \n as a sentinel
*/

#include <stdio.h>

main()
{
	char ch;

	int num_chars = 0;

	while ((ch = getchar()) != '\n')
	{
	putchar(ch);
	num_chars++;
	}

	printf("\nThe number of characters is %d\n", num_chars);
}



