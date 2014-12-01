/* Assignment 2: Question 1 (15 marks) percent_symbol.c */

/* Author: Malachi Griffith
*  Date: Oct. 13 2002
*  Purpose: This program uses nested for statements to print a 
*	    triangular pattern of 'percent' symbols to the screen.
*/

#include <stdio.h>

#define START 6		/* Determines the starting value in the FOR loops. */
#define END 1		/* Determines the ending value in the FOR loops. */

main()
{
	char ch = '%';	/* Variable for the percent symbol so that
			 * escape characters are not needed to display 
	 		 * it in the printf function call.
			*/

	int row;	/* Control variable for the FOR loop. */ 

	int column;	/* Control variable to print '%' symbols in the 
			*  nested FOR loop. */

	int spaces;	/* Control variable to prints spaces in the 
			*  nested FOR loop. */	

	/* Define an "outer" loop 
	*  This loop will start at 6 and iterate until 'row' no longer 
	*  exceeds the value defined by MAXROW (6 times).
	*  This loop will contain two nested loops.  One to print spaces,
	*  and the other to print the '%' symbols. 	
	*/
	
	for (row = START; row >= END; row--)
	{
		printf("\n"); 	/* Brings cursor to a new line before a row */
		
		/* Define the first nested loop.
		*  Prints the character ' ' until 'spaces' is less than
		*  the current value of 'row'.  'spaces' starts at 6 and 
		*  the first value of 'row' is 6, so for the first row, no
		*  spaces are printed.  The 2nd row has 1 space, and so on.
		*/
		for (spaces = START; spaces > row; spaces--)
			printf(" ");
		
		/* Define the second nested loop.
		*  For each loop 'column' is initialized to the value 'row'	
		*  Prints the character 'ch' until 'column' is less than or
		*  equal to the value of 'END'.  For the first row it 
		*  produces 6 '%' symbols, 5 for the second, and so on.
		*/   
		for (column = row; column >= END; column--)
			printf("%c", ch);
	}

	printf("\n\n");

	return(0);
}



