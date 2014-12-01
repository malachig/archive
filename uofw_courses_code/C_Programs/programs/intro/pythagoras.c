/* pythagoras.c*/
/* Author: Malachi Griffith
*  Date: Sept. 20 2002 
*  Purpose: Compute the pythagorean triple (ex. 3, 4, 5 triangle) for any value of m and 
*	    n according to the formula below.
*/

#include <stdio.h>
#include <math.h>

/* Don't forget to compile with cc -lm pythagoras.c */ 

main()
{
 	double side1,
	       side2,
	       m,
	       n,
	       hypotenuse;

	printf ("\n\nEnter the value of m > ");
	scanf (" %lf", &m);
	printf ("Enter the value of n > ");
	scanf (" %lf", &n);

	side1 = (pow(m, 2.0)) - (pow(n, 2.0));
	side2 = 2.0 * m * n;
	hypotenuse = (pow(m, 2.0)) + (pow(n, 2.0));

	printf("\n\nSide 1 = %5.1f\n", side1);
	printf("Side 2 = %5.1f\n", side2); 
	printf("Hypotenuse = %5.1f\n\n", hypotenuse);
	return(0);

}
