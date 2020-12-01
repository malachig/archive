/* Author: Malachi Griffith
*  Date: Oct. 14 2002
*  Purpose: A program for demonstrating the use of \n as a sentinel
*/

#include <stdio.h>

main()
{
	char ch = 'a';

	int num_chars = 0;

	while (ch != '\n')
	{
	ch = getchar();
	putchar(ch);
	num_chars++;
	}

	printf("The number of characters is %d\n", num_chars);
}



