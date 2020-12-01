/* Author: Malachi Griffith
*  Date: Oct. 24 2002
*  Purpose: Sort a data file of prices as an array.
*/

#include <stdio.h>
#define MAX_PRICES 25

/* Function Prototypes */
void readarray(double prices[], int *size);
void sortarray(double prices[]);
void displayarray(double prices[],int size);
void display(double *most_expensive, double *least_expensive);

main()
{
	int size_local;
	double price_list[MAX_PRICES];
	
	readarray(price_list, &size_local);

	displayarray(price_list, size_local); 
	
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
	/*	printf("\nDebug temp is: %f value i=%d", temp, i);
		printf("\nDebug array element is: %f", prices[i]);*/	
		fscanf(input_data, "%lf", &temp); 
		i++;	
	}
	*size = i + 1;	
	fclose(input_data);
}
	
/*
* Function: sortarray
*/




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












