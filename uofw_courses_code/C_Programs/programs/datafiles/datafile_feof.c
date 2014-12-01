/* Author: Malachi Griffith
*  Date: Oct. 1 2002
*  Purpose: Ilustrate the use of the feof() function to specify the end of 
*           the data in a datafile, instead of specifying it as an integer
*	    and then using a for loop.
*/

#include <stdio.h>

main()
{

	int i = 0;

	float data_value;

	FILE * ifptr;

	ifptr = fopen("demo2.dat", "r");
	fscanf(ifptr, " %f", &data_value);

	while(!feof(ifptr))
	{
		i++;
		printf("Data value number %d is: %.2f\n", i, data_value);
		fscanf(ifptr, " %f", &data_value);
	
	}

	printf("There are %d values.\n", i);
	fclose(ifptr);
}
