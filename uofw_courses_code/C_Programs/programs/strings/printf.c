/* Author: Malachi Griffith
*  Date: Nov. 6 2002
*  Purpose: Illustrate the use of te format specifier %s in printf statements.
*  It will start printing from the first position until \0 is encountered.
*/

#include <stdio.h>

#define MAX 30

main()
{
	int j, i;
	char array1[MAX] = "  a b c \0 d e ";

	printf("length of array1 %d\n", strlen(array1));
	printf("%s\n", array1);
}



