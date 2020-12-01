/* Author: Malachi Griffith
*  Date: Oct. 14 2002
*  Purpose: A program to demonstrate some of the collating sequence
*/

#include <stdio.h>

main()
{
	char ch;

	printf("Upper case characters: \n");

	for (ch = 'A'; ch <= 'E'; ch++)
		printf("%3c", ch);
	printf("\n\n");
	printf("appear in positions: \n");

	for (ch = 'A'; ch <= 'E'; ch++)
		printf("%3d", ch);
	printf("\n");

	printf("\n\n%d\n", 'c' < 'A');
}



