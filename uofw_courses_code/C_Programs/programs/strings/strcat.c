/* Author: Malachi Griffith
*  Date: Nov. 6 2002 
*  Purpose: Ilustrate the use of the <string.h> function, strcat().
*  strcat() appends str2 to str1, i.e., str1 is changed and str2 is 
*  unchanged.  Both strings must be null terminated and the result is null
*  terminated.  Str1 must be large enough to accomodate str2 + itself.
*/

#include <stdio.h>
#include <string.h>

#define SIZE1 21
#define SIZE2 11

main()
{
	char str1[SIZE1];
	char str2[SIZE2];

	strcpy(str1, "Hello");
	strcpy(str2, " there!");
	
	strcat(str1, str2);

	printf("%s", str1);
	
	printf("\n");
}



