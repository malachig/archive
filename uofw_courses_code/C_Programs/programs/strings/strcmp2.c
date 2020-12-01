/* Author: Malachi Griffith
*  Date: Nov. 6 2002  
*  Purpose: Illustrates the use of the string compare, strcmp() function.
*  This function compares two strings character by character and returns
*  '0' if they are the same, a positive integer if str1 > str2, and a 
*  negative integer if str1 < str2.
*  This function is in the <string.h> library.
*/

#include <stdio.h>
#include <string.h>

#define SIZE 81

main()
{
	char str1[SIZE];
	char str2[SIZE];

	gets(str1);
	gets(str2);
	printf("The strings are: %d %d\n", strlen(str1), strlen(str2));

	if(!strcmp(str1, str2))
		printf("The strings are equal.\n");
	else
	{ 
		printf("The strings are not equal.\n");
		printf("strcmp gives %d\n", strcmp(str1, str2));
	}
}



