/* Author: Malachi Griffith
*  Date: Oct. 20 2002
*  Purpose: Use a function to take two integer arrays as input and
*  find the absolute difference between values in corresponding locations.
*/

#include <stdio.h>
#include <math.h>

/* Function Prototype */
void absolute_dif(int arr1[], int arr2[], int arr_out[]);

main()
{
	int i;

	int array1[3] = {5, -1, 7};
	int array2[3] = {2, 4, -2};
	int array_out[3];

	absolute_dif(array1, array2, array_out);

	for (i = 0; i <= 2; i++)
		printf("\nThe absolute difference of %d and %d is %d",
		       array1[i], array2[i], array_out[i]);
	
	printf("\n\n");
}

/* Function: absolute_dif */

void
absolute_dif(int arr1[], int arr2[],int arr_out[])
{
	int i;

	for (i = 0; i <= 2; i++)
		arr_out[i] = abs(arr1[i] - arr2[i]);
}
