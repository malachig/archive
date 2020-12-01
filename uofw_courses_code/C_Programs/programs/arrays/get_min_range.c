/* Author: Malachi Griffith
*  Date: Nov. 10 2002
*  Purpose: Returns the subscript of the the smallest value in a portion 
*  	    of an array.  In this case the array contains integer values.
*   	    It has three arguments, an array as input which will not be 
*  	    altered by the function, and a subscript value for the beginning
*	    and ending of the sub-array.
*/

#include <stdio.h>
#define MAX_SIZE 81

/* Function prototype */
int get_min_range(const int data[], int start, int end);


main()
{
	int value;
	int i = 0;
	int array_N[MAX_SIZE] = {0};
	int start, end;
	int min_pos;

	printf("\nAt the prompt please enter the integer data for the array");
	printf("\nHit <enter> after each data point (max of 81 points)");
	printf("\nHit <ctrl> <d> to end. \n");

	while(!feof(stdin))
	{
		scanf("%d", &array_N[i]);
		i++;
	}
	
	printf("\nEnter the starting point for the sub-array > ");
	scanf("%d", &start);
	
	printf("Enter the ending point for the sub-array > ");
	scanf("%d", &end);

	min_pos = get_min_range(array_N, start - 1, end - 1);	
	
	printf("\nThe minimum value in this range is at %d\n", min_pos);
	printf("It's magnitude is %d\n", array_N[min_pos - 1]);
}

int 
get_min_range(const int data[], int start, int end)
{
	int i;
	int small_sub = data[start];
	int min_position = start;

	for (i = start; i <= end; i++)
	{
		if (data[i] < small_sub)
		{
			min_position = i;
			small_sub = data[i];
		}
	}

	return(min_position + 1);
}





