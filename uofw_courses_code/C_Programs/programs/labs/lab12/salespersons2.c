/* salespersons2.c */

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

/* Declare pointer to structure */
typedef sales_record *sales_record_ptr;

/* Function prototypes */
void input(sales_record_ptr data[]);
void sort_averages(sales_record_ptr data[]);
int find_max(sales_record_ptr data[], int start, int end);
void display(sales_record_ptr data[]);
 
main()
{
	sales_record_ptr data[MAX_PEOPLE];
	
	input(data);
	sort_averages(data);
	display(data);

}

/*
*  Function: input()
*/
void 
input(sales_record_ptr data[])
{
	int i, j;	
	sales_record temp;
	int sum = 0;
	FILE *input_data;

	input_data = fopen("salespersons.dat", "r");

	for (i = 0; i < MAX_PEOPLE; i++)
	{
		/* Dynamic Memory Allocation */
		data[i] = (sales_record_ptr)malloc(sizeof(sales_record));	
		fscanf(input_data, "%s", temp.first_name);
		fscanf(input_data, "%s", temp.last_name);

		for (j = 0; j < MAX_YEARS; j++)
		{
			fscanf(input_data, "%d", &temp.sales[j]);
			sum += temp.sales[j];
		}
		temp.average = sum / MAX_YEARS;

		sum = 0;

		*data[i] = temp;
	}
	fclose(input_data);
}

/*
*  Function: sort_averages()
*/
void 
sort_averages(sales_record_ptr data[])
{
	sales_record_ptr temp;
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

/*
* Function: find_max()
*/
int
find_max(sales_record_ptr data[], int start, int end)
{
	int current;
	int i;

	current = start;

	for (i = start + 1; i < end; i++)
		if(data[i]->average > data[current]->average)
			current = i;

	return (current);
}
 
/*
*  Function: display()
*/
void 
display(sales_record_ptr data[])
{
	int i;

	/* Display the average for each salesperson along with their name */
	for (i = 0; i < MAX_PEOPLE; i++)
	{
		printf("\n");
		printf("%s ", data[i]->first_name);
		printf("%s", data[i]->last_name);

		printf("\nThis person's average is %d\n", data[i]->average);
	}
	printf("\n\n");	
}
