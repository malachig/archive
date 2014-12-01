/* temps.c */

/* Author: Malachi Griffith
*  Date: Oct. 26 2002 
*  Purpose: Process a list of temperatures, display each along with its 
*  category (hot, pleasant or cold). As well as print the hottest and coldest
*  temperature.  And finally print the average.
*/

#include <stdio.h>

main()
{
	int temp[] = {55,63,68,74,59,45,41,58,60,67,65,78,82,88,91,92,90,
		      93,87,80,78,79,72,68,61,59};

	int i;
	int highest;
	int lowest;
	int cold = 0, pleasant = 0, hot = 0;

	int N = 26;
	int sum = 0;	
	double average;
	
	highest = temp[0];
	lowest = temp[0];

	printf("\n");

	for(i = 0; i < N; i++)
	{
		if (temp[i] < 60)
		{
			printf("\nDay %d was cold (a temp of %d)",
				 i+1, temp[i]);
			cold++;
		}

		else if (temp[i] <= 84)
		{
			printf("\nDay %d was pleasant (a temp of %d)",	
				i+1, temp[i]);
			pleasant++;
		}

		else if (temp[i] >= 85)
		{
			printf("\nDay %d was hot (a temp of %d)",
				i+1, temp[i]);
			hot++;
		}
		sum += temp[i];

		if (temp[i] > highest)
			highest = temp[i];
		if (temp[i] < lowest)
			lowest = temp[i];
	}
	average = (double) (sum / N);

	printf("\n\nThe hottest day was %d degrees, the coldest %d degrees",
		highest, lowest);
	printf("\nThere were %d cold, %d pleasant, and %d hot days",
		cold, pleasant, hot);
	printf("\nThe average temperature was %.1f\n\n", average);
}



