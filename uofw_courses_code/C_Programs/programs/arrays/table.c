/* Author: Malachi Griffith
*  Date: Oct. 19 2002 
*  Purpose: Store a series of 10 integers in an array, then display them
*  in at table with the percentage that each value is of the total of all
*  values.
*/

#include <stdio.h>

main()
{
	int number[10] = {8,12,18,25,24,30,28,22,23,10};
	double percent[10];
	double sum = 0;
	int i;  /* Control variable for loops */
	char ch = '%';

	printf("\n\nn\t%c of total", ch);
	
	for(i = 1; i <= 10; i++)
		sum += number[i - 1];
		
	for(i = 1; i <= 10; i++)
		percent[i - 1] = (number[i - 1] / sum) * 100;

	for(i = 1; i <= 10; i++)
		printf("\n%d\t%.2f", number[i - 1], percent[i - 1]);
		
	printf("\n\n");

	return(0);
}



