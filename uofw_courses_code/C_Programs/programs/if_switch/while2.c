/* while.c */

/* Author: Malachi Griffith
*  Date: Sept. 20 2002
*  Purpose: Calculate class averages using counter-controlled repetition.
*/

#include <stdio.h>

main()
{
	int counter,
	    class_size;

	float mark,
	      accumulator,
	      average;

	/* Initialise variable */

	counter = 0;
	accumulator = 0.0;

	/* Read input */
	
	printf("\n\nEnter the size of the class > ");
	scanf(" %d", &class_size);

	while (counter < class_size)
		{
		printf("\nPlease enter the mark > ");
		scanf(" %f", &mark);
		accumulator += mark;
		counter++;
		}
	/* Calculate the average */
	
	average = accumulator / (float) class_size;
	printf("\nClass average is %.2f\n\n", average);
}



