/* Author: Malachi Griffith
*  Date: Oct. 24 2002
*  Purpose: Search an array for specified values.
*/

#include <stdio.h>
#define MAX_SIZE 20 

/* Function Prototypes */
void readarray(int values[], int *size);
int searcharray(int values[], int size, int search_value);
void queryarray(int values[], int size);

main()
{
	int size_local;

	int value_list[MAX_SIZE];
	
	readarray(value_list, &size_local);
	
	queryarray(value_list, size_local);
}

/*
* Function readarray
*/
void
readarray(int values[], int *size)
{
	int i = 0;
	double temp;
	
	FILE *input_data;
	input_data = fopen("array.dat", "r");

	fscanf(input_data, "%lf", &temp);

	while(!feof(input_data))
	{
		values[i] = temp;
		fscanf(input_data, "%lf", &temp); 
		i++;	
	}
	*size = i;	
	fclose(input_data);
}
	
/*
*  Function: queryarray
*/
void 
queryarray(int values[], int size)
{
	int target;
	int target_position;

	printf("\nWhich integer value would you like to search for > ");
	scanf("%d", &target);

	target_position = searcharray(values, size, target);

	if (target_position == -1)
		printf("\nThat value is not in the list of data\n\n");
	else
		printf("\nThat value is found at position %d in the array\n\n",
			target_position);	

}

/*
*  Function: searcharray
*/
int 
searcharray(int values[], int size, int search_value)
{
	int i;
	int position = -1;
	
	for(i = 0; i < size; i++)
		if (values[i] == search_value)
			position = i;
	return(position);
}
	
	
	





