/* dowhile.c */

/* Author: Malachi Griffith
*  Date: Sept. 20 2002
*  Purpose: Demonstrates sentinel-controlled repetition using do-while structure
*/

#include <stdio.h>

main()
{
	int counter;
	float 	average,
		mark,
		total_marks;

	/* Initialise variables */
	
	total_marks = 0.0;
	counter = -1;
	mark = 0.0;

	printf("\n\n");
	/* Input data */

	do
	{
		total_marks += mark;
		counter++;
		printf("Please enter mark, or negative number to end: ");
		scanf(" %f", &mark);
	}
	while(!(mark<0.0));

	/* Average data */
	if (counter != 0)
	{
		average = total_marks / counter;
		printf("Class average is %.2f\n\n", average);
	}
	else 
		printf("No marks were entered.\n\n");
}



