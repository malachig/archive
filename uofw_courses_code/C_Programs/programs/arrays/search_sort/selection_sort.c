/* Author: Malachi Griffith
*  Date: Nov. 03 2002  
*  Purpose: Illustrate the basic structure of a selection sort algorithm
*/

#include <stdio.h>
#define SIZE 8

int Find_Min(int[], int, int);

main()
{
	int list[SIZE] = {-3,25,0,14,55,-9,1,43};
	int current;
	int index;
	int temp;

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
		printf("%d, ", list[index]);
	printf("\n\n");
}

int
Find_Min(int list[SIZE], int start, int end)
{
	int current;
	int index;

	current = start;
	for (index = start + 1; index < end; index++)
	  if(list[index] < list[current])
		current = index;
	  return (current);
}


