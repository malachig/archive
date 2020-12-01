/* Author: Malachi Griffith
*  Date: Dec. 15 2002
*  Purpose: Simple use of selection sort
*/

#include <stdio.h>

/* function prototypes */
void CreateArray(int numbers[], int *size);
void PrintArray(int numbers[], int size);
void SortArray(int numbers[], int size);
int FindMin(int numbers[], int start, int end); 

main()
{
	int numbers[20];
	int size_local;

	/* Create Array */
	CreateArray(numbers, &size_local);	
	
	/* Print Array */
	PrintArray(numbers, size_local);

	/* Sort Array */
	SortArray(numbers, size_local);

	/*Print Array */
	PrintArray(numbers, size_local);
}

void 
CreateArray(int numbers[], int *size)
{
	int i = 0;
	FILE *input_data;

	input_data = fopen("numbers.dat", "r");
	
	fscanf(input_data, "%d", &numbers[i]);

	while(!feof(input_data))
	{
		i++;
		fscanf(input_data, "%d", &numbers[i]);
		
	}
	*size = i;
	fclose(input_data);
}

void 
PrintArray(int numbers[], int size)
{
	int i;

	for (i = 0; i < size; i++)
	{
		printf("\n%d", numbers[i]);
	}
	printf("\n");
}



void 
SortArray(int numbers[], int size)
{
	int index;
	int current_min;
	int temp;
	
	for (index = 0; index < size; index++)
	{
		current_min = FindMin(numbers, index, size);

		if (index != current_min)
		{
		  temp = numbers[index];
		  numbers[index] = numbers[current_min];
		  numbers[current_min] = temp;
		}
	}	
}

int 
FindMin(int numbers[], int start, int end)
{
	int i;
	int min_pos;

	min_pos = start;
	
	for(i = start + 1; i < end; i++)
	{
		if (numbers[i] < numbers[min_pos])
			min_pos = i;
	}
	return(min_pos);
}


