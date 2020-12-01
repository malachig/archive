/* Author: Malachi Griffith
*  Date: Dec. 15 2002
*  Purpose: Simple use of selection sort
*/

#include <stdio.h>

/* function prototypes */
void CreateArray(int numbers[], int *size);
void PrintArray(int numbers[], int size);
void SortArray(int numbers[], int size);

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
	int i, j;
	int hold;

	for (i = 1; i < size; i++)
	{
		if (numbers[i] < numbers[i-1]);
		{
		  hold = numbers[i];
		  
		  for (j = i - 1; j >= 0; j--)
		  {
		    if (numbers[j] > hold)
		    	numbers[j+1] = numbers[j];
		    
		    else
			break;
		  }
		  
		  numbers[j+1] = hold;
		}
	}		
}
