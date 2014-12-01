/* Author: Malachi Griffith
*  Date: Oct. 14 2002
*  Purpose: Using EOF to indicate end of data at standard input (keyboard)
*  The key strokes to indicate end of data at the terminal is <ctrl> d.
*/

#include <stdio.h>

main()
{
	char ch;
	printf("Please enter a line of characters\n");
	printf(" followed by pressing enter.\n");
	printf(" enter control d when you want to quit\n");

	while ((ch = getchar()) != EOF)
		printf("%c", ch);
	printf("That is all, folks!\n");
}
