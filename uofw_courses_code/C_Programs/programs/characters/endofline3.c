/* Author: Malachi Griffith
*  Date: Oct. 14 2002
*  Purpose: Use of end of line mark '\n'
*/

#include <stdio.h>

main()
{
	char ch;
	printf("Please enter a line of characters\n");
	printf(" followed by pressing enter.\n");
	while ((ch = getchar()) != '\n')
		printf("%c", ch);
	printf("\n");
	printf("That is all, folks!\n");
}



