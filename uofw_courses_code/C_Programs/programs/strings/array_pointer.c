/* Author: Malachi Griffith
*  Date: Nov. 6 2002
*  Purpose: Array elements are stored in memory in a contiguous fashion.
*  An array name is a constant pointer, containing the address of the first
*  element of the array.  This pointer cannot be changed.  To access all 
*  the elements you can thus make use of a second pointer variable.
*/

#include <stdio.h>
#include <string.h>

#define SIZE 21

main()
{
	char str[SIZE] = "Hello";
	char *str_ptr;

	int i;

	printf("%s\n", str);

	for(i = strlen(str) - 1; i >= 0; i--)
		printf("%c", str[i]);

	printf("\n");

	str_ptr = str + strlen(str) - 1;  /* point to end of the string */

	while (str_ptr >= str)
	{
		putchar(*str_ptr);
		str_ptr--;
	}
	
	printf("\n");
}



