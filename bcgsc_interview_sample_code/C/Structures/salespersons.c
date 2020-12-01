/* salespersons.c */

/* Author: Malachi Griffith
*  Date: Nov. 14 2002 
*  Purpose: Read values describing salepeople from a file and enter the 
*  info into a structure.
*/

#include <stdio.h>
#include <string.h>

#define NAME_SIZE 15
#define MAX_YEARS 5
#define MAX_PEOPLE 4

typedef struct{
	char first_name[NAME_SIZE];
	char last_name[NAME_SIZE];
	int sales[MAX_YEARS];
	} sales_record;


main()
{
	int i, j;	
	sales_record data[MAX_PEOPLE];
	int average[MAX_PEOPLE];
	int sum = 0;
	
	FILE *input_data;
	input_data = fopen("salespersons.dat", "r");

	for (i = 0; i < MAX_PEOPLE; i++)
	{
		fscanf(input_data, "%s", data[i].first_name);
		fscanf(input_data, "%s", data[i].last_name);

		for (j = 0; j < MAX_YEARS; j++)
			fscanf(input_data, "%d", &data[i].sales[j]);
	}
	
	/* Calculate the average yearly sales for each person */	
	for (i = 0; i < MAX_PEOPLE; i++)
	{
		for (j = 0; j < MAX_YEARS; j++)
		{
			sum += data[i].sales[j];
		}
		average[i] = sum / MAX_YEARS;
		sum = 0;
	}


	/* Display the average for each salesperson along with their name */
	for (i = 0; i < MAX_PEOPLE; i++)
	{
		printf("\n");
		printf("%s ", data[i].first_name);
		printf("%s", data[i].last_name);

		printf("\nThis person's average is %d", average[i]);
	}
	printf("\n\n");	



	

	
}



