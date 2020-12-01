/* Author: Malachi Griffith
*  Date: Nov. 6 2002 
*  Purpose: Illustrate the short-comings of the fscanf function when dealing
*  with strings.  The format specifier %s is fine for printf but when used
*  in fscanf it will start inputting at the first non-blank character, and
*  it will stop at the first blank character following the first series of 
*  characters.  (you get only one word).
*/

#include <stdio.h>
#define maximum 10

main()
{
	int j, i;
	char string[maximum];

	FILE *afptr;

	afptr = fopen("data.dat", "r");
	fscanf(afptr, "%s", string);
	
	printf(" string is of length %d\n", strlen(string));
	printf("%s\n", string);
}



