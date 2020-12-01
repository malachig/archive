/* Author: Malachi Griffith
*  Date: Oct. 14 2002
*  Purpose: A program for demonstrating character input and output.
*/

#include <stdio.h>
#define LIMIT 5

main()
{
	char ch;
	int i;

	for (i = 0; i < LIMIT; i++)
	{
		scanf("%c", &ch);
		printf("The character you entered was: %c\n", ch);
	}
}
