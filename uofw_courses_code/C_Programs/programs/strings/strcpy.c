/* Author: Malachi Griffith
*  Date: Nov. 6 2002
*  Purpose: Illustrate the use of the <string.h> function, strcpy().
*  This function copies the contents of str2 into str1.  Str1 must be 
*  large enough to hold str2 or other data may be overwritten.
*/

#include <stdio.h>
#include <string.h>

#define SIZE 81

main()
{
	char str[SIZE];

	strcpy(str, "Hello");

	printf("%s\n", str);
}



