/* Author: Malachi Griffith
*  Date: Nov. 03 2002  
*  Purpose: Illustrate the basic structure of a selection sort algorithm
*/

#include <stdio.h>
#define SIZE 8

int Find_Min(double[], int, int);

main()
{
	double list[SIZE] = {-3.1,25.5,0.6,14.3,55.3,-9.0,1.1,43.1};
	int current;
	int index;
	double temp;

	for (index = 0; index < SIZE; index++)
	{
	   current = Find_Min(list, index, SIZE);
	      if (index != current)
	      {
	      temp = list[index];
	      list[index] = list[current];
	      list[current] = temp;
    	      }
 	}
	printf("\n");
	
	for (index = 0; index < SIZE; index++)
		printf("%.1f, ", list[index]);
	printf("\n\n");
}

int
Find_Min(double list[SIZE], int start, int end)
{
	int current;
	int index;

	current = start;
	for (index = start + 1; index < end; index++)
	  if(list[index] < list[current])
		current = index;
	  return (current);
}


