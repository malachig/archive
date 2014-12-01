/* Author: Malachi Griffith
*  Date: Nov. 03 2002  
*  Purpose: Illustrate the basic structure of a selection sort algorithm
*  This version finds the largest value in a subset and places it in the
*  position n-1, then finds the next largest by scanning an array that
*  ignores the final position and placing it in position n-2, etc.
*  Fills the sorted values from right to left until the beginning of the 
*  array is reached.
*/

#include <stdio.h>
#define SIZE 8

int Find_Max(int[], int, int);

main()
{
	int list[SIZE] = {-3,25,0,14,55,-9,1,43};
	int current;
	int index;
	int temp;

	for (index = SIZE; index >= 0; index--)
	{
	   current = Find_Max(list, index, 0);
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
Find_Max(int list[SIZE], int start, int end)
{
	int current;
	int index;

	current = start;
	for (index = start; index >= end; index--)
	  if(list[index] > list[current])
		current = index;
	  return (current);
}


