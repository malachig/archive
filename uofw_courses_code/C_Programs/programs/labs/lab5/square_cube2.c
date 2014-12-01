/* Assignment 1: Question 2 (10 Marks) */

/* Author: Malachi Griffith
*  Date: Sept. 24 2002 
*  Purpose: To illustrate the use of the while repetition structure.
*           A while loop will be used to calculate the square and cube 
*           of numbers specified by the user.
*  Note: Because this program uses the 'pow' function it must include the
*  math.h library and therefore must be compiled with cc -lm square_cube.c
*/

#include <stdio.h>
#include <math.h>

main()
{
	int number, 	/* Number to be squared and cubed.  Also acts as
			 * the control variable for the loop.         */ 
	    max,	/* Number of loops to perform */		 
      	    square,     /* Squared value of the number.		      */
	    cube;	/* Cubed value of the number.		      */

	/* Ask the user how many squares and cubes to calculate */
	printf("\n\nPlease enter how many numbers you wish to find squares");
	printf("and cubes for > ");
	scanf(" %d", &max);

	/* Initialise the number varible */
	number = 1;

	/* Print the header statement */
	printf("\n\tNumber \tSquare \tCube\n");
	
	/* Use a for loop to find the squares and cubes */
	for(number = 1; number <= max; number++)  
		{
		square = pow(number, 2);
		cube = pow(number, 3);
		printf("\t%d \t%d \t%d\n", number, square, cube);	
		}
	printf("\n");	
	return(0);
}



