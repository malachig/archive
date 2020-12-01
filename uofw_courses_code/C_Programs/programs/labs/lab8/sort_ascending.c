/* Author: Malachi Griffith
*  Date: Oct. 24 2002
*  Purpose: Sort a data file of prices as an array.
*/

#include <stdio.h>
#define MAX_PRICES 25

/* Function Prototypes */
void readarray(double prices[], int *size);
void sortarray(double prices[], int size, double *most_exp, double *least_exp);
double find_min(double prices[], int start, int end);
void displayarray(double prices[],int size);
void display(double most_expensive, double least_expensive);

main()
{
	int size_local;
	double most, least;

	double price_list[MAX_PRICES];
	
	printf("\nPrices Before Sorting of Array:");	
	readarray(price_list, &size_local);
	printf("\n");

	sortarray(price_list, size_local, &most, &least);
	
	printf("\nPrices After Sorting of Array:");
	displayarray(price_list, size_local); 
	printf("\n");

	display(most, least);
	
	printf("\n\n");
}

/*
* Function readarray
*/
void
readarray(double prices[], int *size)
{
	int i = 0;
	double temp;
	
	FILE *input_data;
	input_data = fopen("sort_array.dat", "r");

	fscanf(input_data, "%lf", &temp);

	while(!feof(input_data))
	{
		prices[i] = temp;
		printf("\nThe price: %d is $%.2f", i, prices[i]);	
		fscanf(input_data, "%lf", &temp); 
		i++;	
	}
	*size = i;	
	fclose(input_data);
	printf("\nSize is %d\n", *size);
}
	
/*
* Function: sortarray
*/
void
sortarray(double prices[], int size, double *most_exp, double *least_exp)
{
	int current;
	int index;
	double temp;
	
	*most_exp = prices[0];  /* initialise max and min values */
	*least_exp = prices[0];

	/* consider one value from the array at a time */
	for(index = 0; index <= size - 1; index++)
	{
		/* find current position of minimum value in array */
		/* only search beyond the current position being considered,
		*  ie. the index value. */

		current = find_min(prices, index, size);
		if(index != current)
		
		/* If the min value does not equal the value being considered,
		*  the switch the position with the min value */
		{
		   temp = prices[index];
		   prices[index] = prices[current];
		   prices[current] = temp;
		}
		
		if (prices[index] > *most_exp)
			*most_exp = prices[index];
		if (prices[index] < *least_exp)
			*least_exp = prices[index];
	}
} 


/*
*  Subprogram: find_min 
*/
double
find_min(double prices[], int start, int end)

/* The whole array can be accessed by this function but each time it is
*  called by the sortarray function it is told to start searching from 
*  the index value.  Everything before the index value has already been 
*  sorted.  */

{
	int current;
	int index;

	current = start;  /* Start is index value from the calling function */
			  /* current saves the position of the min value. */

	for (index = start + 1; index < end; index++)
		if(prices[index] < prices[current])
			current = index;
	return (current);
}

/*
* Function: Displayarray
*/
void
displayarray(double prices[], int size)
{
	int i;

	for(i = 0; i < size; i++)
		printf("\nThe price %d: is $%.2f", i, prices[i]);
}

/*
* Function: Display
*/

void 
display(double most_expensive, double least_expensive)
{
	printf("\nThe most expensive item was: $%.2f", most_expensive);
	printf("\nThe least expensive item was: $%.2f\n\n", least_expensive);
}










