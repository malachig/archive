/* Author: Malachi Griffith
*  Date: Nov. 6 2002
*  Purpose: Illustrate the use of the <ctype.h> function toupper().
*  This function converts a string to uppercase.  The <ctype.h> library is
*  the character classification/conversion library. 
*/

#include <stdio.h>
#include <string.h>
#include <ctype.h>

#define SIZE 81

main()
{
	char str[SIZE];
	int i;
	int length;

	strcpy(str, "this is a test");
	length = strlen(str);

	for(i = 0; i < length; i++)
		str[i] = toupper(str[i]);
	printf("%s", str);
	printf("\n");
}



