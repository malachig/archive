/* Author: Malachi Griffith
*  Date: Oct. 3 2002
*  Purpose: read input from the file lab4q2.dat, print the values 
*	    and find their sum.
*/

#include <stdio.h>

main()
{
	int counter = 1;
	int sum = 0;
	int data_value = 0;

	FILE * ifptr;

	ifptr = fopen("lab4q2.dat", "r");
	printf("\n\n");

	while(!feof(ifptr))
	{
	fscanf(ifptr, "%d", &data_value);
	printf("Data value number %d is: %d\n", counter, data_value);
	sum += data_value; 
	counter++;	
	}
	
	printf("\nThe sum of the data values is: %d\n\n", sum);
}



