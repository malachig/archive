/* Author: Malachi Griffith
*  Date: Oct. 8 2002
*  Purpose: Demonstrate the simple usage of data files in programming.
*/

#include <stdio.h>

main()
{
	int i;
	int num_recs;
	
	float data_value;
	FILE * ifptr;	/* The name 'ifptr' is the name given to point to the datafile. */

	ifptr = fopen("demo.dat", "r");  /* The 'r' specifies that it will read from the file */

	fscanf(ifptr, "%d", &num_recs);
	
	for (i = 0; i < num_recs; i++)
	{
		fscanf(ifptr, "%f", &data_value);
		printf("data value number %d is: %.2f\n", i + 1, data_value);
	}

	fclose(ifptr);
}



