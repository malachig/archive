/* Author: Malachi Griffith
*  Date: Nov. 6 2002 
*  Purpose: Uses the exit(1) function from <stdlib.h> to exit from the 
*  program gracefully if there is a problem opening an input file.
*/

#include <stdio.h>
#include <stdlib.h>

#define MAX_CHAR 51

main()
{
	FILE *input_file;
	char file_name[MAX_CHAR];

/*  Get the file name quasi-interactively */

	gets(file_name);

/* Attempt to open the file; generate an error message and 
*  exit gracefully if not possible. */

	if ((input_file = fopen(file_name, "r")) == NULL)
	{
		printf("Error: File %s could not be opened.\n", file_name);
		
		exit(1);
	}
	printf("The file was successfully opened.\n");

/* Attempt to close the file; generate an error message and exit
*  gracefully if not possible */

	if(fclose(input_file) == EOF)
	{
		printf("Error: The file could not be close.\n");
		exit(1);
	}

	printf("The file was successfully closed.\n");
}



