/* Author: Malachi Griffith
*  Date:  Oct. 14 2002
*  Purpose: Processing file for one line of characters using the 
*	    fgetc function.
*/

#include <stdio.h>

main()
{
	char ch;
	
	FILE * infile;
	
	infile = fopen("fgetc.dat", "r");

	while ((ch = fgetc(infile)) != '\n')
		printf("%c", ch);

	printf("\n");

	printf("That is all, folks!\n");
}



