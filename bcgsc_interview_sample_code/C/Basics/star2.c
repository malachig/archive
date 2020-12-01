/* Assignment 1: Question 4 Part B. (10 Marks) star2.c */

/* Author: Malachi Griffith
*  Date: Sept. 28 2002 
*  Purpose: To use a set of nested for statements to display a pattern 
*	    of asterisk characters.
*/

#include <stdio.h>
#define MAXROW 5 	/* Definition of how many rows of stars to print */


main()
{
	int row;	/* Control variable for the following loops */
	int column;	/* Control variable which prints the stars */


	/* Define the 'outer' loop */
	/* Starts at 1 and executes the lines in the braces following
	*  until 'row' exceed the MAXROW value of 5. */

	for (row = 1; row <= MAXROW; row++)
	{	
		printf("\n\t");
		
		/* Define the nested loop. */
		/* Prints stars according to the 'column' control variable.
		*  When the column exceeds the number of the current row 
		*  the loop breaks and control is returned to the outer 	
		*  loop.  Therefore, for row 1 it produces 1 star, 2 stars
		*  for row 2, and so on.
		*/ 	
	
		for (column = 1; column <= row; column++)
		printf("*");
	}

	printf("\n\n");
	
	return(0);
}



