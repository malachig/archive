/* Author: Malachi Griffith
*  Date: Nov. 6 2002
*  Purpose: Give functions for safely opening and closing multiple
*  input and output files.  Problems will result in error messages and 
*  graceful exit from the program.
*/

#include <stdio.h>
#include <stdlib.h>

#define ITER 20

/*Function Prototypes */

FILE * SafeFopen (char *, char *);
void SafeFclose (FILE *);

main()
{
	FILE *out1_file;
	FILE *out2_file;

	int i;
	int num;

/* Open the two output files */

	out1_file = SafeFopen("m1.dat", "w");
	out2_file = SafeFopen("m2.dat", "w");

/*  Use a random number generator to create numbers for output.
*  Write odd number to file m1.dat, and the even numbers to m2.dat.
*/

	for (i = 1; i < ITER; i++)
	{
		num = rand();
		if (num %2)
			fprintf(out1_file, "%d\n", num);
		else
			fprintf(out2_file, "%d\n", num);
	}

/* Close the files */
	
	SafeFclose(out1_file);
	SafeFclose(out2_file);
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



