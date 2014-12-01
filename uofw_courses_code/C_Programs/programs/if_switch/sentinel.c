/* sentinel.c */

/* Author: Malachi Griffith
*  Date: Sept. 20 2002
*  Purpose: Calculate class averages using sentinel-controlled repitition.
*/

#include <stdio.h>

main()
{
	int counter;

	float average,
	      mark,
	      total_marks;

	/* Initialise variables */

	total_marks = 0.0;
	counter = 0;

	/* Input data */

	printf("\n\nPlease enter mark, or negative number to end > ");
	scanf(" %f", &mark);
	
	while(!(mark < 0.0))
		{
		total_marks += mark;
		counter++;
		printf("Please enter mark, or negative number to end > ");
		scanf(" %f", &mark);
		}
	/* Average data and print output */

	if(counter != 0)
		{
		average = total_marks / (float) counter;
		printf("Class average is %.2f\n\n", average);
		}
	else
		printf("No marks were entered.\n\n");

}



