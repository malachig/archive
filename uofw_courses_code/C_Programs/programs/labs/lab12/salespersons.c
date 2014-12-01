/* Author: Malachi Griffith
*  Date: Nov. 14 2002 
*  Purpose: Read values describing salepeople from a file and enter the 
*  info into a structure. NOW SORT IN DESCENDING ORDER BY AVERAGE SALES
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
	int average;
	} sales_record;

/* Function prototypes */
void sort_averages(sales_record data[]);
int find_max(sales_record data[], int start, int end);

main()
{
	int i, j;	
	sales_record data[MAX_PEOPLE];
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
		data[i].average = sum / MAX_YEARS;
		sum = 0;
	}


	sort_averages(data);

	/* Display the average for each salesperson along with their name */
	for (i = 0; i < MAX_PEOPLE; i++)
	{
		printf("\n");
		printf("%s ", data[i].first_name);
		printf("%s", data[i].last_name);

		printf("\nThis person's average is %d", data[i].average);
	}
	printf("\n\n");	
}

void 
sort_averages(sales_record data[])
{
	sales_record temp;
	int index;	
	int current;
	int min_position;


	for (index = 0; index < MAX_PEOPLE - 1; index++)
	{
		current = find_max(data, index, MAX_PEOPLE);

		if(index != current)
		{
			temp = data[index];
			data[index] = data[current];
			data[current] = temp;
		}
	}
}	

int
find_max(sales_record data[], int start, int end)
{
	int current;
	int i;

	current = start;

	for (i = start + 1; i < end; i++)
		if(data[i].average > data[current].average)
			current = i;

	return (current);
}
 








