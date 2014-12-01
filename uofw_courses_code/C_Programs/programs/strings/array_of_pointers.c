/* Author: Malachi Griffith
*  Date: Nov. 18 2002 
*  Purpose: Orders a list of strings according to string length- shortest 
*  to longest.  The functions used receives an input/output argument that
*  is an array of pointers.
*/

#include <stdio.h>
#include <string.h>

#define MAX_LENGTH 81
#define MAX_STRINGS 20

/* Function Prototype */
void sort_by_size(char *sorted[], int size);
int find_shortest (char *sorted[], int start, int end);

int
main()
{
	char strings[MAX_STRINGS][MAX_LENGTH];

	char *sorted[MAX_STRINGS];

	int num_strings;  /* actual number of strings */
	int i = 0;
	int size;

	/* Get series of strings from user */
	printf("\nEnter each string followed by enter\n");
	
	while(!feof(stdin))
	{
		gets(strings[i]);
		i++;
	}
	size = i;

	/* Fill array of pointers */
	for (i = 0; i <= size; i++)
		sorted[i] = strings[i];
 	
	sort_by_size(sorted, size);  
	
	printf("\n");
	
	for (i = 0; i <= size; i++)
		printf("\n%s", sorted[i]);

	printf("\n");
	
	return(0);
}

/*
*  Function: sort_by_size()
*  Sorts strings by their length
*/
void
sort_by_size(char *sorted[], int size)
{
	int i;
	int index_of_min;

	char *temp;
	
	for (i = 0; i < size; i++)
	{
		index_of_min = find_shortest(sorted, i, size);
	
		if (index_of_min != i)
		{
			/* Sort addresses of 'sorted' */
			temp = sorted[index_of_min];
			sorted[index_of_min] = sorted[i];
			sorted[i] = temp;
		}
	}
}

/*
*  Function: find_shortest()
*/
int
find_shortest (char *sorted[], int start, int end)
{
	int shortest;
	int i;

	shortest = start;

	for (i = start + 1; i <= end; i++)
		if ( (strlen(sorted[i])) < (strlen(sorted[shortest])) )
			shortest = i;

	return(shortest);
}



