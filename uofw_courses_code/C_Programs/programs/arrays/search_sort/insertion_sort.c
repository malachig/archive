/* Author: Malachi Griffith
*  Date: Nov. 2 2002
*  Purpose: To illustrate the basic structure of an inserion sort algorithm.
*/

#include <stdio.h>
#define SIZE 8

main()
{
	int i, j;
	int hold;

	int data[SIZE];

	scanf("%d", &data[0]);

	for (i = 1; i < SIZE; i++)
	{
		scanf("%d", &data[i]);

		if (data[i] < data [i-1])
		{
		   hold = data[i];

		   for(j = i - 1; j >= 0; j--)
		   {
		   	if(data[j] > hold)
				data[j + 1] = data[j];
			else
				break;
	       	   }

		   data[j + 1] = hold;

		   printf("Inserting data point [%d] in position %d \n",
		       i, (j + 1));
		}	
	}
	printf("\n\n");

	for(i = 0; i < SIZE; i++)
		printf("Data[%d] = %d\n", i, data[i]);

}



