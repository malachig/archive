/* Author: Malachi Griffith
*  Date: Nov. 6 2002
*  Purpose: Print the string entered by keyboard in reverse.
*/

#include <stdio.h>
#include <string.h>

#define SIZE 81

main()
{
	char str[SIZE];

	int i;

	gets(str);
	for(i = strlen(str) - 1; i >= 0; i--)
		printf("%c", str[i]);
	printf("\n");
}



