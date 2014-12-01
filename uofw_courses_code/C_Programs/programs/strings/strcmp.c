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

#define SIZE 21

main()
{
	char str1[SIZE] = "marigold";
	char str2[SIZE] = "tulip";

	printf("The value of 'm' = %d\n", 'm');

	printf("The value of 't' = %d\n", 't');

	printf("m - t = %d\n", ('m' - 't'));

	printf("strcmp gives: %d\n", strcmp(str1, str2));

}



