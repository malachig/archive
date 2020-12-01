/* Author: Malachi Griffith
*  Date: Nov. 6 2002 
*  Purpose: Illustrate the use of the <string.h> function strlen().
*  This function returns the length of the string specifed.
*/

#include <stdio.h>
#include <string.h>

#define SIZE 81

main()
{
	char str[SIZE];

	gets(str);

	printf("%d\n", strlen(str));
}



