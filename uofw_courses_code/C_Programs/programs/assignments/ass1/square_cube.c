/* Assignment 1: Question 2 (10 Marks)  square_cube.c */

/* Author: Malachi Griffith
*  Date: Sept. 24 2002 
*  Purpose: To illustrate the use of the while repetition structure.
*           A while loop will be used to calculate the square and cube 
*           of the numbers: 1,2,3, and 4.
*  Note: Because this program uses the 'pow' function it must include
*  the math.h library and therefore must be compiled with,
*  cc -lm square_cube.c
*/

#include <stdio.h>
#include <math.h>

main()
{
	int number, 	/* Number to be squared and cubed.  Also acts as
			 * the control variable for the loop.         */ 
	    square,     /* Squared value of the number.		      */
	    cube;	/* Cubed value of the number.		      */

	/* Initialise the number variable */
	number = 1;

	/* Print the header statement */
	printf("\n\tNumber \tSquare \tCube\n");
	
	/* Use a while loop to find the squares and cubes */
	while (number <= 4)
		{
		square = pow(number, 2);
		cube = pow(number, 3);
		printf("\t%d \t%d \t%d\n", number, square, cube);	
		number++;
		}
	printf("\n");	

	return(0);
}



