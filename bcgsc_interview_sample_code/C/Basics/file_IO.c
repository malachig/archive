/* Author: Malachi Griffith
*  Date: Nov. 6 2002
*  Purpose: Give functions for safely opening and closing multiple
*  input and output files.  Problems will result in error messages and 
*  graceful exit from the program.
*/

#include <stdio.h>
#include <stdlib.h>

/*Function Prototypes */

FILE * SafeFopen (char *, char *);
void SafeFclose (FILE *);

main()
{
	FILE *in1_file;
	FILE *in2_file;

	int i;
	int num;

/* Open the input files */

	in1_file = SafeFopen("m1.dat", "r");
	in2_file = SafeFopen("m2.dat", "r");

/* Read the numbers just written to these files by the previous program,
*  safe_files.c.  Find the remainder resulting from dividing each number
*  by 10, and write it to stdout. */ 

	while(!feof(in1_file))
	{
		fscanf(in1_file, "%d\n", &num);
		num %= 10;
		printf("%d\n", num);
	}

	printf(" ----------------\n");

	while (!feof(in2_file))
	{
		fscanf(in2_file, "%d\n", &num);
		num %= 10;
		printf("%d\n", num);
	}

/* Close the files */
	
	SafeFclose(in1_file);
	SafeFclose(in2_file);
}

FILE 
* SafeFopen(char *f_name, char *mode)
{
	FILE *file_ptr;	

	if ((file_ptr = fopen(f_name, mode)) == NULL)
	{
		fprintf(stderr, "ERROR: File %s could not be opened.\n", 
			f_name);
		exit(1);
	}

	return (file_ptr);
}

void
SafeFclose (FILE *file_ptr)
{
	if (fclose(file_ptr) == EOF)
	{
	fprintf(stderr, "ERROR: File cannot be closed.\n");
		exit(1);
	}
	
	return;
}



