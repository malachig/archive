/* Author: Malachi Griffith
*  Date: Oct. 20 2002 
*  Purpose: Use a function to receive elements from an array and 
*  calculate the product of them.
*/

#include <stdio.h>

/* Function prototype */
void multiply (int array1[], int num_elements, double *product);


main()
{
	int values[10] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
	int num;
	double result_product;

	printf("\nHow many array values would you like to multiply > ");
	scanf("%d", &num);

	multiply(values, num, &result_product);

	printf("\nThe product is %f\n", result_product);
}

/* Function: multiply */
void
multiply (int array1[], int num_elements, double *product)
{
	int i;
	*product = 1;

	for (i = 1; i <= num_elements; i++)
		*product *= (double) array1[i-1];
}
 

