/* Author: Malachi Griffith
*  Date:  Nov. 2 2002 
*  Purpose: Illustrate the basic structure of the bubble sort algorithm
*/

#include <stdio.h>
#define SIZE 8

main()
{
	int list[] = {-3,25,0,14,55,-9,1,43};
	int limit;
	int exchange;
	int index;
	int temp;
	int flag;

	index = 0;
	limit = SIZE - 1;
	flag = 1;
	exchange = 0;

	while(flag)
	{
		flag = 0;
		for(index = 0; index < limit; ++index)
		{
			if (list[index] > list[index + 1])
			{
				temp = list[index];
				list[index] = list[index + 1];
				list[index + 1] = temp;
				flag = 1;
				exchange = index;
			}
		}
	limit = exchange;
	}
	printf("\n\n");

	for (index = 0; index < SIZE; ++index)
		printf("%d, ", list[index]);
	
	printf("\n\n");

}



