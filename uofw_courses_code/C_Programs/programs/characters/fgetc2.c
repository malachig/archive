/* Author: Malachi Griffith
*  Date: Oct. 14 2002 
*  Purpose: Processing a text file until end of data
*/

#include <stdio.h>

FILE * infile;

main()
{
	char ch;
	infile = fopen("fgetc2.dat", "r");
	ch = fgetc(infile);
	while (!feof(infile))
	{
		printf("%c", ch);
		ch = fgetc(infile);
	}

	printf("That is all, folks!\n");
}
