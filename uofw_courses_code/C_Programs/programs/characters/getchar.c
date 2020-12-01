/* Author: Malachi Griffith
*  Date: Oct. 14 2002
*  Purpose: Illustrate use of the getchar command.
*/

#include <stdio.h>

main()
{
	int ch;
	printf("Please enter a character > ");
	ch = getchar();
	
	putchar(ch);
	putchar('\n');
}
