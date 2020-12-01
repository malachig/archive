/* Author: Malachi Griffith
*  Date: Oct. 20 2002
*  Purpose:  Use a function to take an input array argument and find 
*  the absolute value of each element and then print them out as a table
*/

#include <stdio.h>
#include <math.h>

/* Function Prototype */
void absolute(double array1[], int num_elements);

main()
{
	double data;
	double values[25];
	int i = 0;

	printf("\nEnter up to 25 numbers, enter 0 when finished\n");
	
	for (scanf(" %lf", &data);
	     data != 99;
	     scanf(" %lf", &data))
	{
	values[i - 1] = data;
	i++;
	}
	absolute (values, i);
}

/* Function absolute */

void
absolute (double array1[], int num_elements)
{
	int j;
	double abs_value;
	
	printf("\n\tx\t|x|");
	for(j = 0; j < num_elements; j++)
	{
		abs_value = fabs(array1[j-1]);
		printf("\n\t%.1f\t%.1f", array1[j-1], abs_value);
	}	

	printf("\n\n");
}
 


