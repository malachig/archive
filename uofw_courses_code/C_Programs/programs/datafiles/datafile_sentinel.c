/* Author: Malachi Griffith
*  Date: Oct. 1 2002
*  Purpose:  Illustrate the use of a sentinel value to determine
*	     when a program should stop reading data from a file which 
*	     has an unknown number of data in it.
*/

#include <stdio.h>

main()
{
	int counter;
	
	float average;
	float mark;
	float total_marks;
	FILE * ifptr;
	
	/* Initialise the variables */

	total_marks = 0.0;
	counter = 0;
	ifptr = fopen("demo3.dat", "r");

	fscanf(ifptr, "%f", &mark);

	while((mark >= 0.0))
	{
		total_marks += mark;
		counter++;
		fscanf(ifptr, "%f", &mark);
	}

	/* Average data and print output */
	
	if(counter != 0)
	{
		average = total_marks / (float)counter;
		printf("Class average is %.2f\n", average);
		printf("Class has %d students\n", counter);
	}
	else
		printf("No marks were entered.\n");
	
	fclose(ifptr);
}
